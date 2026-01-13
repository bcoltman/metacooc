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
    
    def __init__(self, cap: int = 0, dtype=np.int32):
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
    si, ei = indptr[i], indptr[i + 1]
    sj, ej = indptr[j], indptr[j + 1]
    ni = ei - si
    nj = ej - sj
    if ni == 0 and nj == 0:
        return
        
    row_i = indices[si:ei]
    row_j = indices[sj:ej]
    
    ws.ensure(ni + nj, dtype=indices.dtype)
    
    pi = pj = 0
    n_inter = n_ai = n_aj = 0
    while pi < ni and pj < nj:
        vi = row_i[pi]
        vj = row_j[pj]
        if vi == vj:
            ws.inter[n_inter] = vi
            n_inter += 1
            pi += 1
            pj += 1
        elif vi < vj:
            ws.ai[n_ai] = vi
            n_ai += 1
            pi += 1
        else:
            ws.aj[n_aj] = vj
            n_aj += 1
            pj += 1
            
    if pi < ni:
        tail = row_i[pi:ni]
        ws.ai[n_ai:n_ai + tail.size] = tail
        n_ai += tail.size
    if pj < nj:
        tail = row_j[pj:nj]
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
    nonempty_rows: np.ndarray,
    rng: np.random.Generator,
    n_steps: int,
    ws: _CurveballWorkspace,
) -> None:
    if n_steps <= 0 or nonempty_rows.size < 2:
        return
    n = nonempty_rows.size
    for _ in range(n_steps):
        a = rng.integers(0, n)
        b = rng.integers(0, n - 1)
        if b >= a:
            b += 1
        i = int(nonempty_rows[a])
        j = int(nonempty_rows[b])
        _curveball_trade(indptr, indices, i, j, rng, ws)


def curveball_markov(
    X_csr: sp.csr_matrix,
    n_reps: int,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
) -> Iterable[sp.csr_matrix]:
    """
    FF: curveball on prepared CSR. Yields CSR matrices with dtype int8.
    Preserves row and column sums exactly.
    """
    X0 = X_csr
    rng = _rng_from_state(random_state)
    
    indptr = X0.indptr
    indices = X0.indices
    
    row_nnz = indptr[1:] - indptr[:-1]
    nonempty_rows = np.flatnonzero(row_nnz > 0).astype(np.int64, copy=False)
    
    if burn_in_steps is None:
        burn_in_steps = max(1000, 10 * int(nonempty_rows.size))
    if steps_per_rep is None:
        steps_per_rep = max(int(nonempty_rows.size), 10)
        
    ws = _CurveballWorkspace(cap=0, dtype=indices.dtype)
    _curveball_steps(indptr, indices, nonempty_rows, rng, burn_in_steps, ws)
    
    for rep in range(n_reps):
        if rep > 0:
            _curveball_steps(indptr, indices, nonempty_rows, rng, steps_per_rep, ws)
            
        Y = X0.copy()
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
    burn_in_steps: Optional[int] = None,   # FF only
    steps_per_rep: Optional[int] = None,   # FF only
    sort_indices: bool = False,            # speed knob
) -> Iterable[sp.csr_matrix]:
    """
    Unified interface for null model generation.
    
    Models supported in this reduced module:
      - FF: fixed row + fixed col totals (curveball Markov chain)
      - FE: fixed row totals; columns equiprobable
      - EF: fixed col totals; rows equiprobable
      - EE: equiprobable cells (preserves total fill exactly)
      
    Notes
    -----
    - All samplers yield CSR matrices with dtype=int8 and data=ones.
    - 'preserve_fill' is always ON by design in this reduced module.
    - By default, does NOT sort indices on outputs (often fastest).
      Enable sort_indices=True only if downstream ops benefit.
    """
    if n_reps <= 0:
        return iter(())
        
    model = str(model).upper()
    
    if model in ("EF",):
        Xp = prepare_presence_matrix(X, fmt="csc", copy=copy)
        return ef_equiprob_rows_fixed_cols(Xp, n_reps, random_state=random_state, sort_indices=sort_indices)
        
    Xp = prepare_presence_matrix(X, fmt="csr", copy=copy)
    
    if model == "FF":
        return curveball_markov(
            Xp,
            n_reps=n_reps,
            burn_in_steps=burn_in_steps,
            steps_per_rep=steps_per_rep,
            random_state=random_state,
            sort_indices=sort_indices,
        )
    if model == "FE":
        return fe_fixed_rows_equiprob_cols(Xp, n_reps, random_state=random_state, sort_indices=sort_indices)
    if model == "EE":
        return ee_equiprobable(Xp, n_reps, random_state=random_state, sort_indices=sort_indices)
        
    raise ValueError(f"Unknown or unsupported null model: {model}")




def _col_counts_for_rows_csr(X_csr: sp.csr_matrix, rows: np.ndarray, n_cols: int) -> np.ndarray:
    indptr = X_csr.indptr
    indices = X_csr.indices
    
    nnz = 0
    for r in rows:
        nnz += (indptr[r + 1] - indptr[r])
    if nnz == 0:
        return np.zeros(n_cols, dtype=np.int32)
        
    buf = np.empty(nnz, dtype=indices.dtype)
    pos = 0
    for r in rows:
        s = indptr[r]; e = indptr[r + 1]
        k = e - s
        if k:
            buf[pos:pos + k] = indices[s:e]
            pos += k
            
    return np.bincount(buf, minlength=n_cols).astype(np.int32, copy=False)





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
_G_term_rows = None
_G_nonterm_rows = None
_G_subset_idx = None
_G_n_cols = None
_G_N_T = None

_G_iA = None
_G_iB = None

_G_progress_q = None
_G_progress_every = None

_G_chunk_rows = None
_G_structure_do_nodf = None


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
    global _G_X, _G_model, _G_sort, _G_burn, _G_steps, _G_stat_fn, _G_obs
    global _G_term_rows, _G_nonterm_rows, _G_subset_idx, _G_n_cols, _G_N_T
    global _G_iA, _G_iB, _G_metric_fns
    global _G_chunk_rows, _G_structure_do_nodf
    
    _G_X = X
    _G_model = model
    _G_sort = sort_indices
    _G_burn = burn_in_steps
    _G_steps = steps_per_rep
    _G_stat_fn = stat_fn
    _G_obs = obs
    
    # pull optional items from init_kwargs
    _G_term_rows = init_kwargs.get("term_rows", None)
    _G_nonterm_rows = init_kwargs.get("nonterm_rows", None)
    _G_subset_idx = init_kwargs.get("subset_idx", None)
    _G_n_cols = init_kwargs.get("n_cols", None)
    _G_N_T = init_kwargs.get("N_T", None)
    
    _G_iA = init_kwargs.get("iA", None)
    _G_iB = init_kwargs.get("iB", None)
    
    _G_metric_fns = init_kwargs.get("metric_fns", None)
    
    _G_chunk_rows = init_kwargs.get("chunk_rows", None)
    _G_structure_do_nodf = init_kwargs.get("structure_do_nodf", True)



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
):
    global _G_progress_q, _G_progress_every
    _G_progress_q = progress_q
    _G_progress_every = int(progress_every) if progress_every is not None else 0
    
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


def _worker_run(task):
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
    
    it = null_matrices(
        _G_X, model=_G_model, n_reps=int(n_reps_local),
        random_state=int(seed), copy=False,
        burn_in_steps=_G_burn, steps_per_rep=_G_steps,
        sort_indices=_G_sort,
    )
    
    processed = 0
    every = int(_G_progress_every or 0)
    
    for X_null in it:
        acc["n_done"] += 1
        try:
            X_null = X_null.tocsr()
            v = _G_stat_fn(X_null)
            if v is None:
                acc["n_err"] += 1
            else:
                v = np.asarray(v, dtype=float)
                if v.shape != obs.shape or not np.all(np.isfinite(v)):
                    acc["n_err"] += 1
                else:
                    acc["sum"] += v
                    acc["sum2"] += v * v
                    acc["ge"] += (v >= obs)
                    acc["n_ok"] += 1
        except Exception:
            acc["n_err"] += 1
            
        processed += 1
        if every > 0 and _G_progress_q is not None and (processed % every) == 0:
            _G_progress_q.put(every)
            
    if _G_progress_q is not None:
        if every > 0:
            rem = processed % every
            if rem:
                _G_progress_q.put(rem)
        else:
            _G_progress_q.put(processed)
            
    return acc


def stat_fn_association_jaccard(X_null_full: sp.csr_matrix) -> np.ndarray:
    # Uses worker globals set in _worker_init(...)
    counts_T = _col_counts_for_rows_csr(X_null_full, _G_term_rows, _G_n_cols)
    counts_notT = _col_counts_for_rows_csr(X_null_full, _G_nonterm_rows, _G_n_cols)
    
    a_null = counts_T[_G_subset_idx].astype(np.float64, copy=False)
    b_null = counts_notT[_G_subset_idx].astype(np.float64, copy=False)
    
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(
            a_null,
            (b_null + _G_N_T),
            out=np.zeros_like(a_null, dtype=float),
            where=(b_null + _G_N_T) > 0,
        )


def stat_fn_cooccurrence_jaccard(X_null_full: sp.csr_matrix) -> np.ndarray:
    # Uses worker globals set in _worker_init(...)
    X_null_sub = X_null_full[:, _G_subset_idx].tocsr()
    
    totals = np.asarray(X_null_sub.sum(axis=0)).ravel().astype(np.float64, copy=False)
    
    co = (X_null_sub.T @ X_null_sub).tocsr()
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

def stat_fn_structure_metrics(X_null: sp.csr_matrix) -> np.ndarray | None:
    """
    Returns [c_score, mean_jaccard, nodf] for a null matrix, or None if invalid.

    Uses:
      - _G_chunk_rows (int)
      - _G_structure_do_nodf (bool)
    """
    # Import lazily to avoid circular imports at module import time.
    # These functions must be importable from a stable module path.
    try:
        from metacooc.structure import compute_c_score, mean_jaccard_dot, compute_nodf_streamed
    except Exception:
        # If this triggers circulars, move these metric functions into e.g. metacooc.metrics
        # and import from there instead.
        return None
        
    chunk_rows = int(_G_chunk_rows) if _G_chunk_rows is not None else 50_000
    do_nodf = bool(_G_structure_do_nodf)
    
    try:
        v0 = float(compute_c_score(X_null, chunk_rows=chunk_rows))
        v1 = float(mean_jaccard_dot(X_null, chunk_rows=chunk_rows))
        
        if do_nodf:
            v2 = float(compute_nodf_streamed(X_null, chunk_rows=chunk_rows))
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
    X_prepared: sp.spmatrix,
    model: str,
    n_reps: int,
    obs: np.ndarray,
    stat_fn: Callable[[sp.csr_matrix], np.ndarray],
    random_state: Optional[int] = None,
    n_workers: Optional[int] = None,
    chunksize: int = 1,  # unused in this async/poll pattern; keep for API stability
    sort_indices: bool = False,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    mp_start: str = "fork",
    progress_every: int = 1,  # tqdm granularity per worker
    **init_kwargs,
) -> dict:
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
        n_workers = 4  # or: max(1, (mp.cpu_count() or 1) - 1)
    n_workers = max(1, int(n_workers))
    
    # workers == chains
    reps_per = _split_reps(n_reps, n_workers)
    
    # deterministic per-worker seeds
    seeds = _seed_seq_spawn(random_state, n_workers)
    
    # only launch workers that actually have work
    tasks = [(reps_per[i], seeds[i]) for i in range(n_workers) if reps_per[i] > 0]
    
    n_target = int(sum(r for r, _ in tasks))
    if n_target != n_reps:
        raise RuntimeError(
            f"Requested n_reps={n_reps}, but scheduled n_target={n_target}. "
            f"n_workers={n_workers}, reps_per={reps_per}, tasks={tasks[:5]}..."
        )
        
    ctx = mp.get_context(mp_start)
    progress_q = ctx.Queue()
    
    acc0 = {
        "sum": np.zeros_like(obs, dtype=float),
        "sum2": np.zeros_like(obs, dtype=float),
        "ge": np.zeros_like(obs, dtype=np.int64),
        "n_ok": 0,
        "n_err": 0,
        "n_done": 0,
        "n_target": n_target,
    }
    
    with ctx.Pool(
        processes=min(n_workers, len(tasks)),
        initializer=_worker_init_wrap,
        initargs=(
            X_prepared,
            str(model).upper(),
            bool(sort_indices),
            burn_in_steps,
            steps_per_rep,
            stat_fn,
            obs,
            dict(init_kwargs),
            progress_q,
            int(progress_every),
        ),
    ) as pool:
        async_results = [pool.apply_async(_worker_run, (t,)) for t in tasks]
        
        # tqdm should track what we actually scheduled
        with tqdm(
            total=n_target,
            desc=f"Null ({str(model).upper()}) - each process updates every {progress_every} nulls",
            dynamic_ncols=True,
        ) as pbar:
            import queue as _queue
            import time
            
            # IMPORTANT FIX:
            # Loop until all async_results have been collected.
            # Do NOT use `done < len(async_results)` because async_results shrinks.
            while async_results:
                # ---- drain progress queue ----
                drained_any = False
                try:
                    inc = progress_q.get(timeout=0.1)
                except _queue.Empty:
                    inc = None
                    
                if inc is not None:
                    drained_any = True
                    pbar.update(int(inc))
                    while True:
                        try:
                            inc2 = progress_q.get_nowait()
                        except _queue.Empty:
                            break
                        else:
                            pbar.update(int(inc2))
                            
                # ---- collect finished workers (merge once) ----
                still_pending = []
                any_finished = False
                for r in async_results:
                    if r.ready():
                        any_finished = True
                        acc_local = r.get()
                        acc0 = _merge_accumulators(acc0, acc_local)
                    else:
                        still_pending.append(r)
                async_results = still_pending
                
                # ---- brief yield if nothing happened ----
                if (not drained_any) and (not any_finished):
                    time.sleep(0.05)
                    
    out = _finalise_accumulator(obs, acc0)
    
    # Safety: confirm we merged all scheduled reps
    if int(out.get("n_done", 0)) != int(out.get("n_target", 0)):
        raise RuntimeError(
            f"Internal error: merged n_done={out.get('n_done')} but expected n_target={out.get('n_target')}."
        )
        
    return out