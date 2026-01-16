# null_models.py
from __future__ import annotations

from typing import Callable, Dict, Iterable, Optional, Tuple, Literal

from tqdm import tqdm
import numpy as np
import scipy.sparse as sp
import multiprocessing as mp

SIMMODEL = Literal["FF", "FE", "EF", "EE"]

_OUT_DTYPE = np.int8

def _best_mp_start() -> str:
    """
    Cross-platform start-method chooser:
      - fork if available (Linux; sometimes macOS if explicitly enabled)
      - otherwise spawn (Windows/macOS default)
    """
    methods = mp.get_all_start_methods()
    return "fork" if "fork" in methods else "spawn"

# -----------------------------------------------------------------------------
# Preparation (done once in null_matrices)
# -----------------------------------------------------------------------------

def prepare_presence_matrix(
    X: sp.spmatrix,
    *,
    fmt: Literal["csr", "csc"] = "csr",
    copy: bool = False,
) -> sp.spmatrix:
    """
    One-time normalisation:
      - convert to CSR/CSC
      - eliminate explicit zeros
      - sum duplicates
      - sort indices (canonical)
    """
    if fmt == "csr":
        Xp = X.tocsr(copy=copy)
    elif fmt == "csc":
        Xp = X.tocsc(copy=copy)
    else:
        raise ValueError("fmt must be 'csr' or 'csc'")
    Xp.eliminate_zeros()
    Xp.sum_duplicates()
    Xp.sort_indices()
    return Xp


def _rng_from_state(random_state: Optional[int | np.random.Generator]) -> np.random.Generator:
    return random_state if isinstance(random_state, np.random.Generator) else np.random.default_rng(random_state)


# -----------------------------------------------------------------------------
# FF (Curveball) — prepared CSR in/out (always CSR)
# -----------------------------------------------------------------------------

class _CurveballWorkspace:
    __slots__ = ("inter", "ai", "aj", "pool")
    
    def __init__(self, cap: int = 0, dtype=np.int64):
        self.inter = np.empty(cap, dtype=dtype)
        self.ai = np.empty(cap, dtype=dtype)
        self.aj = np.empty(cap, dtype=dtype)
        self.pool = np.empty(cap, dtype=dtype)
        
    def ensure(self, cap: int, dtype):
        if self.inter.size < cap:
            new_cap = int(max(cap, self.inter.size * 1.5 + 8))
            self.inter = np.empty(new_cap, dtype=dtype)
            self.ai = np.empty(new_cap, dtype=dtype)
            self.aj = np.empty(new_cap, dtype=dtype)
            self.pool = np.empty(new_cap, dtype=dtype)


def _merge_fill(a: np.ndarray, b: np.ndarray, out: np.ndarray) -> None:
    i = j = k = 0
    na, nb = a.size, b.size
    while i < na and j < nb:
        av, bv = a[i], b[j]
        if av <= bv:
            out[k] = av
            i += 1
        else:
            out[k] = bv
            j += 1
        k += 1
    if i < na:
        out[k:k + (na - i)] = a[i:na]
    elif j < nb:
        out[k:k + (nb - j)] = b[j:nb]

def _curveball_trade(
    indptr: np.ndarray,
    indices: np.ndarray,
    i: int,
    j: int,
    rng: np.random.Generator,
    ws: _CurveballWorkspace,
) -> None:
    """
    Curveball trade between two entities i and j, represented by:
      - indptr: offsets into indices
      - indices[indptr[k]:indptr[k+1]]: sorted items in entity k

    Works identically for:
      - CSR rows (entities=rows, items=columns)
      - CSC cols (entities=cols, items=rows)
    """
    si, ei = int(indptr[i]), int(indptr[i + 1])
    sj, ej = int(indptr[j]), int(indptr[j + 1])
    ni = ei - si
    nj = ej - sj
    if ni == 0 and nj == 0:
        return

    a = indices[si:ei]
    b = indices[sj:ej]

    ws.ensure(ni + nj, dtype=indices.dtype)

    pi = pj = 0
    n_inter = n_ai = n_aj = 0

    # a and b must be sorted (prepare_presence_matrix() ensures this)
    while pi < ni and pj < nj:
        va = a[pi]
        vb = b[pj]
        if va == vb:
            ws.inter[n_inter] = va
            n_inter += 1
            pi += 1
            pj += 1
        elif va < vb:
            ws.ai[n_ai] = va
            n_ai += 1
            pi += 1
        else:
            ws.aj[n_aj] = vb
            n_aj += 1
            pj += 1

    if pi < ni:
        tail = a[pi:ni]
        ws.ai[n_ai:n_ai + tail.size] = tail
        n_ai += tail.size

    if pj < nj:
        tail = b[pj:nj]
        ws.aj[n_aj:n_aj + tail.size] = tail
        n_aj += tail.size

    if n_ai == 0 and n_aj == 0:
        return

    ws.pool[:n_ai] = ws.ai[:n_ai]
    ws.pool[n_ai:n_ai + n_aj] = ws.aj[:n_aj]
    pool = ws.pool[:n_ai + n_aj]
    rng.shuffle(pool)

    new_ai = pool[:n_ai]
    new_aj = pool[n_ai:]
    new_ai.sort()
    new_aj.sort()

    inter = ws.inter[:n_inter]
    _merge_fill(inter, new_ai, indices[si:ei])
    _merge_fill(inter, new_aj, indices[sj:ej])

def _curveball_steps(
    indptr: np.ndarray,
    indices: np.ndarray,
    nonempty_entities: np.ndarray,
    rng: np.random.Generator,
    n_steps: int,
    ws: _CurveballWorkspace,
) -> None:
    if n_steps <= 0 or nonempty_entities.size < 2:
        return
    n = int(nonempty_entities.size)
    for _ in range(int(n_steps)):
        a = int(rng.integers(0, n))
        b = int(rng.integers(0, n - 1))
        if b >= a:
            b += 1
        i = int(nonempty_entities[a])
        j = int(nonempty_entities[b])
        _curveball_trade(indptr, indices, i, j, rng, ws)


def curveball_markov(
    X: sp.spmatrix,
    *,
    n_reps: int,
    burn_in_steps: int,
    steps_per_rep: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    Curveball Markov chain on a prepared sparse matrix.
    
    Behaviour depends solely on X:
      - CSR  → swaps columns between row pairs (species-based)
      - CSC  → swaps rows between column pairs (site-based)
      
    X must already be prepared (CSR or CSC, sorted indices).
    """
    rng = _rng_from_state(random_state)
    
    indptr = X.indptr
    indices = X.indices
    
    nnz_per_entity = indptr[1:] - indptr[:-1]
    nonempty = np.flatnonzero(nnz_per_entity > 0).astype(np.int64, copy=False)
    
    ws = _CurveballWorkspace(dtype=indices.dtype)
    
    _curveball_steps(indptr, indices, nonempty, rng, burn_in_steps, ws)
    
    for rep in range(n_reps):
        if rep > 0:
            _curveball_steps(indptr, indices, nonempty, rng, steps_per_rep, ws)
            
        Y = X.tocsr(copy=True)
        Y.data = np.ones(Y.indices.size, dtype=_OUT_DTYPE)
        if sort_indices:
            Y.sort_indices()
        yield Y


# -----------------------------------------------------------------------------
# FE / EF / EE — direct samplers (all preserve fill; all yield CSR int8)
# -----------------------------------------------------------------------------

def fe_fixed_rows_equiprob_cols(
    X_csr: sp.csr_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    FE: fixed row totals, columns equiprobable.
    """
    n_rows, n_cols = X_csr.shape
    row_deg = (X_csr.indptr[1:] - X_csr.indptr[:-1]).astype(np.int64, copy=False)
    N = int(row_deg.sum())
    
    indptr = np.empty(n_rows + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(row_deg, out=indptr[1:])
    
    data = np.ones(N, dtype=_OUT_DTYPE)
    rng = _rng_from_state(random_state)
    
    for _ in range(n_reps):
        indices = np.empty(N, dtype=np.int64)
        pos = 0
        for i in range(n_rows):
            k = int(row_deg[i])
            if k:
                indices[pos:pos + k] = rng.choice(n_cols, size=k, replace=False)
                pos += k
        Y = sp.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
        if sort_indices:
            Y.sort_indices()
        yield Y


def ef_equiprob_rows_fixed_cols(
    X_csc: sp.csc_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    EF: fixed column totals, rows equiprobable.
    Construct in CSC then convert to CSR per replicate.
    """
    n_rows, n_cols = X_csc.shape
    col_deg = (X_csc.indptr[1:] - X_csc.indptr[:-1]).astype(np.int64, copy=False)
    N = int(col_deg.sum())
    
    indptr = np.empty(n_cols + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(col_deg, out=indptr[1:])
    
    data = np.ones(N, dtype=_OUT_DTYPE)
    rng = _rng_from_state(random_state)
    
    for _ in range(n_reps):
        indices = np.empty(N, dtype=np.int64)
        for j in range(n_cols):
            k = int(col_deg[j])
            if k:
                s = int(indptr[j]); e = int(indptr[j + 1])
                indices[s:e] = rng.choice(n_rows, size=k, replace=False)
        Y = sp.csc_matrix((data, indices, indptr), shape=(n_rows, n_cols)).tocsr()
        if sort_indices:
            Y.sort_indices()
        yield Y


def ee_equiprobable(
    X_csr: sp.csr_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    EE:
    - preserves total fill exactly (N = nnz(X))
    - samples N distinct cells uniformly without replacement
    """
    n_rows, n_cols = X_csr.shape
    N = int(X_csr.nnz)
    n_cells = int(n_rows) * int(n_cols)
    
    rng = _rng_from_state(random_state)
    
    if N == 0 or n_cells == 0:
        for _ in range(n_reps):
            yield sp.csr_matrix((n_rows, n_cols), dtype=_OUT_DTYPE)
        return
        
    if N > n_cells:
        raise ValueError("nnz exceeds total number of cells; invalid input matrix.")
        
    for _ in range(n_reps):
        lin = rng.choice(n_cells, size=N, replace=False)
        
        rows = (lin // n_cols).astype(np.int64, copy=False)
        cols = (lin % n_cols).astype(np.int64, copy=False)
        
        # Build CSR without COO->CSR:
        order = np.lexsort((cols, rows))
        rows = rows[order]
        cols = cols[order]
        
        row_counts = np.bincount(rows, minlength=n_rows).astype(np.int64, copy=False)
        indptr = np.empty(n_rows + 1, dtype=np.int64)
        indptr[0] = 0
        np.cumsum(row_counts, out=indptr[1:])
        
        data = np.ones(N, dtype=_OUT_DTYPE)
        Y = sp.csr_matrix((data, cols, indptr), shape=(n_rows, n_cols))
        
        if sort_indices:
            # already sorted by construction, but keep for consistency
            Y.sort_indices()
            
        yield Y


# -----------------------------------------------------------------------------
# Public entry point — prepares once; ALL outputs are CSR int8
# -----------------------------------------------------------------------------

def null_matrices(
    X: sp.spmatrix,
    model: SIMMODEL,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    copy: bool = False,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    Unified interface for null model generation.
    
    Models supported in this reduced module:
      - FF: fixed row + fixed col totals (curveball Markov chain)
      - FE: fixed row totals; columns equiprobable
      - EF: fixed col totals; rows equiprobable
      - EE: equiprobable cells (preserves total fill exactly)
    """
    
    if n_reps <= 0:
        return iter(())

    model = str(model).upper()

    if model == "FF":
        # X is already prepared and oriented by the parent
        return curveball_markov(
            X,
            n_reps=n_reps,
            burn_in_steps=burn_in_steps,
            steps_per_rep=steps_per_rep,
            random_state=random_state,
            sort_indices=sort_indices,
        )

    if model == "EF":
        return ef_equiprob_rows_fixed_cols(X, n_reps, random_state=random_state, sort_indices=sort_indices)

    if model == "FE":
        return fe_fixed_rows_equiprob_cols(X, n_reps, random_state=random_state, sort_indices=sort_indices)

    if model == "EE":
        return ee_equiprobable(X, n_reps, random_state=random_state, sort_indices=sort_indices)

    raise ValueError(f"Unknown or unsupported null model: {model}")



def _row_counts_for_cols_csr(X_csr: sp.csr_matrix, col_mask: np.ndarray, out_len: int) -> np.ndarray:
    # Reuse preallocated _G_sel_buf
    indptr = X_csr.indptr
    sel = col_mask[X_csr.indices].astype(np.int64, copy=False)  # view, no new array if bool->int64 forces copy; unavoidable cast
    # pack into preallocated buffer
    buf = _G_sel_buf
    buf[:sel.size] = sel
    buf[sel.size] = 0  # padding for reduceat on empty rows
    counts = np.add.reduceat(buf, indptr[:-1])
    if counts.size != out_len:
        counts = counts[:out_len]
    # Fix empty rows: reduceat on padding already yields 0, so usually no-op
    return counts

def _seed_seq_spawn(master_seed: Optional[int], n: int) -> list[int]:
    """
    Deterministic per-worker seeds.
    Uses SeedSequence so streams are independent even if adjacent seeds.
    """
    ss = np.random.SeedSequence(0 if master_seed is None else int(master_seed))
    return [int(s.generate_state(1, dtype=np.uint32)[0]) for s in ss.spawn(n)]


def _split_reps(n_reps: int, n_workers: int) -> list[int]:
    n_workers = max(1, int(n_workers))
    base = n_reps // n_workers
    rem = n_reps % n_workers
    return [base + (1 if i < rem else 0) for i in range(n_workers)]

def _chunk_reps(n_reps: int, block: int) -> list[int]:
    """
    Split n_reps into chunks of size <= block.
    """
    out = []
    while n_reps > 0:
        k = min(block, n_reps)
        out.append(k)
        n_reps -= k
    return out

# ---- core accumulator merge ----

def _merge_accumulators(a: dict, b: dict) -> dict:
    # Both are dicts with:
    #  sum: ndarray, sum2: ndarray, ge: ndarray, n_ok: int, n_err: int
    a["sum"] += b["sum"]
    a["sum2"] += b["sum2"]
    a["ge"] += b["ge"]
    a["n_ok"] += int(b["n_ok"])
    a["n_err"] += int(b["n_err"])
    a["n_done"] += int(b["n_done"])
    return a


def _finalise_accumulator(obs: np.ndarray, acc: dict) -> dict:
    """
    Turn sum/sum2/ge into mean/sd/ses/p_emp.
    Matches your (k+1)/(n+1) convention.
    """
    n_ok = int(acc["n_ok"])
    out = {
        "mean": np.full_like(obs, np.nan, dtype=float),
        "sd": np.full_like(obs, np.nan, dtype=float),
        "ses": np.full_like(obs, np.nan, dtype=float),
        "p_emp": np.full_like(obs, np.nan, dtype=float),
        "n_ok": n_ok,
        "n_err": int(acc["n_err"]),
        "n_done": int(acc["n_done"]),
        "n_target": int(acc["n_target"]),
    }
    if n_ok <= 0:
        return out
        
    mean = acc["sum"] / n_ok
    var = acc["sum2"] / n_ok - mean * mean
    var[var < 0] = 0.0
    sd = np.sqrt(var)
    
    out["mean"] = mean
    out["sd"] = sd
    
    with np.errstate(divide="ignore", invalid="ignore"):
        out["ses"] = np.divide(
            obs - mean,
            sd,
            out=np.full_like(obs, np.nan, dtype=float),
            where=sd > 0,
        )
        
    out["p_emp"] = (acc["ge"] + 1.0) / (n_ok + 1.0)
    return out



# -----------------------------------------------------------------------------
# Parallel “reduce over null replicates” for vector-valued statistics
# -----------------------------------------------------------------------------

# Worker globals to avoid pickling big objects repeatedly
_G_X = None
_G_model = None
_G_sort = None
_G_burn = None
_G_steps = None
_G_stat_fn = None
_G_obs = None

# --- worker globals for stat_fns ---
_G_term_cols = None
_G_nonterm_cols = None
_G_subset_idx = None
_G_n_rows = None
_G_N_T = None

_G_iA = None
_G_iB = None

_G_progress_q = None
_G_progress_every = None

_G_chunk_rows = None
_G_structure_do_nodf = None

# ---- FF chain state (per worker) ----
_G_ff_rng = None
_G_ff_ws = None
_G_ff_nonempty = None
_G_ff_indptr = None
_G_ff_indices = None

# ---- burn-in progress plumbing (parent -> workers) ----
_G_burn_q = None
_G_burn_every = 0
_G_seed_q = None


def _ff_chain_init_and_burn() -> None:
    """
    FF-only: initialise a persistent Curveball chain *once per worker* and perform
    burn-in with coarse progress updates to the parent via _G_burn_q.

    Assumes:
      - _G_X is already prepared CSR or CSC with sorted indices.
      - _G_burn and _G_steps are set.
      - _G_seed_q provides one seed per worker.
    """
    global _G_ff_rng, _G_ff_ws, _G_ff_nonempty, _G_ff_indptr, _G_ff_indices

    # Get a deterministic per-worker seed from the seed queue
    seed = None
    if _G_seed_q is not None:
        seed = int(_G_seed_q.get())
    _G_ff_rng = np.random.default_rng(seed)

    X = _G_X
    _G_ff_indptr = X.indptr
    _G_ff_indices = X.indices

    nnz_per_entity = _G_ff_indptr[1:] - _G_ff_indptr[:-1]
    _G_ff_nonempty = np.flatnonzero(nnz_per_entity > 0).astype(np.int64, copy=False)

    _G_ff_ws = _CurveballWorkspace(dtype=_G_ff_indices.dtype)

    burn = int(_G_burn or 0)
    if burn <= 0 or _G_ff_nonempty.size < 2:
        # signal "done" for this worker
        if _G_burn_q is not None:
            _G_burn_q.put(("DONE", 1))
        return

    every = int(_G_burn_every or 0)

    if every <= 0 or _G_burn_q is None:
        # No progress reporting; just do it
        _curveball_steps(_G_ff_indptr, _G_ff_indices, _G_ff_nonempty, _G_ff_rng, burn, _G_ff_ws)
        if _G_burn_q is not None:
            _G_burn_q.put(("DONE", 1))
        return

    # Coarse-grained progress reporting: do burn-in in blocks
    done = 0
    while done < burn:
        step = min(every, burn - done)
        _curveball_steps(_G_ff_indptr, _G_ff_indices, _G_ff_nonempty, _G_ff_rng, step, _G_ff_ws)
        done += step
        _G_burn_q.put(int(step))

    _G_burn_q.put(("DONE", 1))



def _worker_init(
    X: sp.spmatrix,
    model: str,
    sort_indices: bool,
    burn_in_steps: Optional[int],
    steps_per_rep: Optional[int],
    stat_fn,          # callable at module scope
    obs: np.ndarray,
    **init_kwargs,
):
    """
    Fast worker initialiser:
      - caches matrix, parameters, stat_fn, and obs
      - builds reusable boolean masks for term / non-term columns
      - preallocates a selection buffer of size (nnz + 1) for reduceat
    """
    global _G_X, _G_model, _G_sort, _G_burn, _G_steps, _G_stat_fn, _G_obs
    global _G_term_cols, _G_nonterm_cols, _G_subset_idx, _G_n_rows, _G_N_T
    global _G_iA, _G_iB, _G_metric_fns
    global _G_chunk_rows, _G_structure_do_nodf

    # new caches for accelerated _row_counts_for_cols_csr
    global _G_col_mask_T, _G_col_mask_notT, _G_sel_buf, _G_sel_buf_len

    _G_X = X
    _G_model = model
    _G_sort = sort_indices
    _G_burn = burn_in_steps
    _G_steps = steps_per_rep
    _G_stat_fn = stat_fn
    _G_obs = obs

    # pull optional items from init_kwargs
    _G_term_cols     = init_kwargs.get("term_cols", None)
    _G_nonterm_cols  = init_kwargs.get("nonterm_cols", None)
    _G_subset_idx    = init_kwargs.get("subset_idx", None)
    _G_n_rows        = init_kwargs.get("n_rows", None)
    _G_N_T           = init_kwargs.get("N_T", None)

    _G_iA            = init_kwargs.get("iA", None)
    _G_iB            = init_kwargs.get("iB", None)

    _G_metric_fns    = init_kwargs.get("metric_fns", None)

    _G_chunk_rows    = init_kwargs.get("chunk_rows", None)
    _G_structure_do_nodf = init_kwargs.get("structure_do_nodf", True)

    # ---- precompute term / non-term masks when available ----
    _G_col_mask_T = None
    _G_col_mask_notT = None
    if _G_term_cols is not None or _G_nonterm_cols is not None:
        n_cols = int(X.shape[1])
        if _G_term_cols is not None:
            _G_col_mask_T = np.zeros(n_cols, dtype=np.bool_)
            _G_col_mask_T[np.asarray(_G_term_cols, dtype=np.int64)] = True
        if _G_nonterm_cols is not None:
            _G_col_mask_notT = np.zeros(n_cols, dtype=np.bool_)
            _G_col_mask_notT[np.asarray(_G_nonterm_cols, dtype=np.int64)] = True

    # ---- preallocate selection buffer for reduceat ----
    # size = nnz + 1 (padding sentinel for empty-row safety)
    nnz = int(X.nnz)
    _G_sel_buf_len = nnz + 1
    _G_sel_buf = np.empty(_G_sel_buf_len, dtype=np.int64)


def _worker_init_wrap(
    X: sp.spmatrix,
    model: str,
    sort_indices: bool,
    burn_in_steps: Optional[int],
    steps_per_rep: Optional[int],
    stat_fn,
    obs: np.ndarray,
    init_kwargs: dict,
    progress_q,
    progress_every: int,
    burn_q,
    burn_every: int,
    seed_q,
):
    global _G_progress_q, _G_progress_every
    global _G_burn_q, _G_burn_every, _G_seed_q
    
    # keep existing (currently unused)
    _G_progress_q = progress_q
    _G_progress_every = int(progress_every) if progress_every is not None else 0
    
    # burn-in plumbing
    _G_burn_q = burn_q
    _G_burn_every = int(burn_every) if burn_every is not None else 0
    _G_seed_q = seed_q
    
    _worker_init(
        X,
        model,
        sort_indices,
        burn_in_steps,
        steps_per_rep,
        stat_fn,
        obs,
        **(init_kwargs or {}),
    )
    
    # FF: initialise persistent chain and burn-in once per worker
    if str(model).upper() == "FF":
        _ff_chain_init_and_burn()

def _worker_run_chunk(task):
    """
    Run a chunk of null replicates within a worker.

    task: (chunk_reps: int, chunk_seed: int)

    Notes:
      - For FF, chunk_seed is ignored because the worker runs a persistent chain
        seeded once in _worker_init_wrap, and burn-in was done once per worker.
      - For FE/EF/EE, chunk_seed is used (direct sampling, no chain).
    """
    n_reps_local, seed = task
    obs = _G_obs

    acc = {
        "sum": np.zeros_like(obs, dtype=float),
        "sum2": np.zeros_like(obs, dtype=float),
        "ge": np.zeros_like(obs, dtype=np.int64),
        "n_ok": 0,
        "n_err": 0,
        "n_done": 0,
    }

    # ---- FF path: persistent chain per worker ----
    if str(_G_model).upper() == "FF":
        reps = int(n_reps_local)
        steps = int(_G_steps or 0)

        for _ in range(reps):
            # advance chain between nulls
            if steps > 0:
                _curveball_steps(_G_ff_indptr, _G_ff_indices, _G_ff_nonempty, _G_ff_rng, steps, _G_ff_ws)

            # materialise current state as CSR for downstream stats
            Y = _G_X.tocsr(copy=True)
            Y.data = np.ones(Y.indices.size, dtype=_OUT_DTYPE)
            if _G_sort:
                Y.sort_indices()

            acc["n_done"] += 1
            try:
                v = _G_stat_fn(Y)
                if v is None:
                    acc["n_err"] += 1
                    continue

                v = np.asarray(v, dtype=float)
                if v.shape != obs.shape or not np.all(np.isfinite(v)):
                    acc["n_err"] += 1
                    continue

                acc["sum"] += v
                acc["sum2"] += v * v
                acc["ge"] += (v >= obs)
                acc["n_ok"] += 1

            except Exception:
                acc["n_err"] += 1

        return acc, reps

    # ---- non-FF path (FE/EF/EE): direct sampler per chunk seed ----
    it = null_matrices(
        _G_X,
        model=_G_model,
        n_reps=int(n_reps_local),
        random_state=int(seed),
        copy=False,
        burn_in_steps=_G_burn,
        steps_per_rep=_G_steps,
        sort_indices=_G_sort,
    )

    for X_null in it:
        acc["n_done"] += 1
        try:
            v = _G_stat_fn(X_null)
            if v is None:
                acc["n_err"] += 1
                continue

            v = np.asarray(v, dtype=float)
            if v.shape != obs.shape or not np.all(np.isfinite(v)):
                acc["n_err"] += 1
                continue

            acc["sum"] += v
            acc["sum2"] += v * v
            acc["ge"] += (v >= obs)
            acc["n_ok"] += 1

        except Exception:
            acc["n_err"] += 1

    return acc, int(n_reps_local)


def stat_fn_association_jaccard(X: sp.csr_matrix) -> np.ndarray:
    counts_T    = _row_counts_for_cols_csr(X, _G_col_mask_T, _G_n_rows).astype(np.float64, copy=False)
    counts_notT = _row_counts_for_cols_csr(X, _G_col_mask_notT, _G_n_rows).astype(np.float64, copy=False)

    a_null = counts_T[_G_subset_idx]
    b_null = counts_notT[_G_subset_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(a_null, (b_null + _G_N_T),
                         out=np.zeros_like(a_null, dtype=float),
                         where=(b_null + _G_N_T) > 0)


def stat_fn_cooccurrence_jaccard(X: sp.csr_matrix) -> np.ndarray:
    # New convention: X is taxa × samples, and _G_subset_idx indexes taxa (rows).
    X_sub = X[_G_subset_idx, :].tocsr()
    
    totals = np.asarray(X_sub.sum(axis=1)).ravel().astype(np.float64, copy=False)
    
    co = (X_sub @ X_sub.T).tocsr()
    inter = co[_G_iA, _G_iB].A1.astype(np.float64, copy=False)
    
    A = totals[_G_iA]
    B = totals[_G_iB]
    union = A + B - inter
    
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(
            inter,
            union,
            out=np.zeros_like(inter, dtype=float),
            where=union > 0,
        )

def stat_fn_structure_metrics(X: sp.csr_matrix) -> np.ndarray | None:
    """
    Returns [c_score, mean_jaccard, nodf] for a null matrix (taxa × samples), or None if invalid.
    
    Uses worker globals:
      - _G_chunk_rows (int)
      - _G_structure_do_nodf (bool)
      
    """
    from metacooc.structure import compute_c_score, mean_jaccard_dot, compute_nodf_streamed
    
    chunk_rows = int(_G_chunk_rows) if _G_chunk_rows is not None else 50_000
    do_nodf = bool(_G_structure_do_nodf)
    
    try:
        v0 = float(compute_c_score(X, chunk_rows=chunk_rows))
        v1 = float(mean_jaccard_dot(X, chunk_rows=chunk_rows))
        
        if do_nodf:
            v2 = float(compute_nodf_streamed(X, chunk_rows=chunk_rows))
        else:
            v2 = np.nan
            
        if not (np.isfinite(v0) and np.isfinite(v1) and (np.isfinite(v2) or not do_nodf)):
            return None
            
        return np.array([v0, v1, v2], dtype=float)
        
    except MemoryError:
        return None
    except Exception:
        return None


def parallel_null_reduce_vector(
    *,
    X: sp.spmatrix,
    model: str,
    n_reps: int,
    obs: np.ndarray,
    stat_fn: Callable[[sp.csr_matrix], np.ndarray],
    random_state: Optional[int] = None,
    n_workers: Optional[int] = None,
    chunksize: int = 1,  # Pool.imap_unordered chunksize
    sort_indices: bool = False,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    mp_start: str = "fork",
    progress_every: int = 1,  # interpreted as preferred chunk size
    **init_kwargs,
) -> dict:
    """
    Parallel reduction over null replicates using chunked tasks and imap_unordered.

    FF behaviour:
      - Matrix is prepared once in the parent, oriented CSR/CSC based on shape.
      - Each worker runs an independent Markov chain:
          - seeded once
          - burn-in once (with parent burn-in tqdm)
          - then emits nulls by stepping steps_per_rep between samples
      - Chunk seeds are ignored for FF (chain-based, not chunk-based).

    FE/EF/EE behaviour:
      - Matrix is prepared once in parent (CSR/CSC as needed).
      - Chunk seeds are used for deterministic direct sampling.
    """
    n_reps = int(n_reps)
    if n_reps <= 0:
        return {
            "mean": np.full_like(obs, np.nan, dtype=float),
            "sd": np.full_like(obs, np.nan, dtype=float),
            "ses": np.full_like(obs, np.nan, dtype=float),
            "p_emp": np.full_like(obs, np.nan, dtype=float),
            "n_ok": 0,
            "n_err": 0,
            "n_done": 0,
            "n_target": 0,
        }

    if n_workers is None:
        n_workers = 4
    n_workers = max(1, int(n_workers))

    model_u = str(model).upper()

    # ---- chunk plan ----
    chunk_n = int(progress_every) if progress_every and progress_every > 0 else 50
    chunk_n = max(1, chunk_n)

    chunk_sizes = _chunk_reps(n_reps, block=chunk_n)

    # deterministic per-chunk seeds (used for non-FF; ignored for FF)
    ss = np.random.SeedSequence(0 if random_state is None else int(random_state))
    child_seqs = ss.spawn(len(chunk_sizes))
    chunk_seeds = [int(cs.generate_state(1, dtype=np.uint32)[0]) for cs in child_seqs]

    tasks = list(zip(chunk_sizes, chunk_seeds))
    n_target = int(sum(chunk_sizes))

    ctx = mp.get_context(mp_start)

    acc0 = {
        "sum": np.zeros_like(obs, dtype=float),
        "sum2": np.zeros_like(obs, dtype=float),
        "ge": np.zeros_like(obs, dtype=np.int64),
        "n_ok": 0,
        "n_err": 0,
        "n_done": 0,
        "n_target": n_target,
    }

    # ---- prepare matrix once in the parent ----
    if model_u == "FF":
        n_rows, n_cols = map(int, X.shape)

        if n_rows <= n_cols:
            X = prepare_presence_matrix(X, fmt="csr", copy=False)
            mode_msg = "INFO: Curveball on rows (species); swapping samples between species"
        else:
            X = prepare_presence_matrix(X, fmt="csc", copy=False)
            mode_msg = "INFO: Curveball on columns (samples); swapping species between samples"

        print(mode_msg)

        basis = max(n_rows, n_cols)
        if burn_in_steps is None:
            burn_in_steps = max(1000, 5 * basis)
        if steps_per_rep is None:
            steps_per_rep = max(basis, 10)

        print(f"INFO: Performing {int(burn_in_steps)} burn-in shuffles per chain before first statistics.")

    elif model_u == "EF":
        X = prepare_presence_matrix(X, fmt="csc", copy=False)

    elif model_u in ("FE", "EE"):
        X = prepare_presence_matrix(X, fmt="csr", copy=False)

    else:
        raise ValueError(f"Unknown or unsupported null model: {model_u}")

    from tqdm import tqdm as _tqdm
    import queue as _queue

    n_procs = min(n_workers, len(tasks))

    # ---- FF burn-in progress setup ----
    burn_q = None
    seed_q = None
    burn_every = 0

    if model_u == "FF":
        burn_q = ctx.Queue()
        seed_q = ctx.Queue()

        # one seed per worker (chain seed), derived deterministically from master seed
        worker_seeds = _seed_seq_spawn(random_state, n_procs)
        for s in worker_seeds:
            seed_q.put(int(s))

        # coarse progress emission every N trades per worker
        # tune as needed; larger = less IPC, smaller = smoother bar
        burn_every = 2000

    # ---- launch pool ----
    with ctx.Pool(
        processes=n_procs,
        initializer=_worker_init_wrap,
        initargs=(
            X,
            model_u,
            bool(sort_indices),
            burn_in_steps,
            steps_per_rep,
            stat_fn,
            obs,
            dict(init_kwargs),
            None,  # progress_q unused
            0,     # progress_every unused
            burn_q,
            burn_every,
            seed_q,
        ),
    ) as pool:

        # ---- show burn-in tqdm (FF only) ----
        if model_u == "FF":
            burn_steps = int(burn_in_steps or 0)
            if burn_steps > 0:
                total = burn_steps * n_procs
                done_workers = 0

                with _tqdm(
                    total=total,
                    desc=f"Burn-in (FF) - {n_procs} chains",
                    dynamic_ncols=True,
                ) as pbar_burn:
                    while done_workers < n_procs:
                        try:
                            msg = burn_q.get(timeout=0.2)
                        except _queue.Empty:
                            continue

                        if isinstance(msg, tuple) and len(msg) == 2 and msg[0] == "DONE":
                            done_workers += int(msg[1])
                        else:
                            pbar_burn.update(int(msg))

        # ---- main null tqdm ----
        desc = f"Null ({model_u}) - chunks of {chunk_n}"
        with _tqdm(total=n_target, desc=desc, dynamic_ncols=True) as pbar:
            for acc_local, chunk_k in pool.imap_unordered(
                _worker_run_chunk,
                tasks,
                chunksize=max(1, int(chunksize)),
            ):
                acc0 = _merge_accumulators(acc0, acc_local)
                pbar.update(int(chunk_k))

    out = _finalise_accumulator(obs, acc0)

    if int(out.get("n_done", 0)) != int(out.get("n_target", 0)):
        raise RuntimeError(
            f"Internal error: merged n_done={out.get('n_done')} but expected n_target={out.get('n_target')}."
        )

    return out

