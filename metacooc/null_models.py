# null_models.py
from __future__ import annotations

from typing import Callable, Iterable, Optional, Literal, Protocol

import multiprocessing as mp
import numpy as np
import scipy.sparse as sp


SIMMODEL = Literal["FF", "FE", "EF", "EE"]
_OUT_DTYPE = np.int8


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def _best_mp_start() -> str:
    """
    Cross-platform start-method chooser:
      - fork if available (Linux; sometimes macOS if explicitly enabled)
      - otherwise spawn (Windows/macOS default)
    """
    methods = mp.get_all_start_methods()
    return "fork" if "fork" in methods else "spawn"
    
def _rng_from_state(random_state: Optional[int | np.random.Generator]) -> np.random.Generator:
    """
    Return a NumPy Generator from an integer seed or an existing Generator.
    """
    return random_state if isinstance(random_state, np.random.Generator) else np.random.default_rng(random_state)


def _record_debug_event(debug_callback, **event) -> None:
    """
    Send a sampler decision to an optional benchmark/debug sink.
    """
    if debug_callback is None:
        return
    payload = dict(event)
    if callable(debug_callback):
        debug_callback(payload)
    elif hasattr(debug_callback, "append"):
        debug_callback.append(payload)
def prepare_presence_matrix(
    X: sp.spmatrix,
    *,
    fmt: Literal["csr", "csc"] = "csr",
    copy: bool = False,
) -> sp.spmatrix:
    """
    Convert a sparse matrix to CSR/CSC and canonicalise internal structure.

    Operations:
      - convert to CSR/CSC
      - eliminate explicit zeros
      - sum duplicates
      - sort indices
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


# -----------------------------------------------------------------------------
# FF (Curveball
# -----------------------------------------------------------------------------

class _CurveballWorkspace:
    __slots__ = ("inter", "ai", "aj", "pool")

    def __init__(self, cap: int = 0, dtype=np.int64):
        self.inter = np.empty(cap, dtype=dtype)
        self.ai = np.empty(cap, dtype=dtype)
        self.aj = np.empty(cap, dtype=dtype)
        self.pool = np.empty(cap, dtype=dtype)

    def ensure(self, cap: int, dtype) -> None:
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
    Perform a Curveball trade between two entities i and j.

    indptr/indices encode sorted item sets for each entity:
      - items of entity k are indices[indptr[k]:indptr[k+1]]

    Works for both:
      - CSR (entities=rows, items=columns)
      - CSC (entities=columns, items=rows)
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
    """
    Apply n_steps random Curveball trades among non-empty entities.
    """
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


# -----------------------------------------------------------------------------
# FE / EF / EE direct generators
# -----------------------------------------------------------------------------

def _sparse_indices_dtype(n_items: int) -> np.dtype:
    """
    Smallest SciPy-compatible index dtype for item indices.
    """
    return np.int32 if int(n_items) <= np.iinfo(np.int32).max else np.int64


def _sparse_indptr_dtype(nnz: int) -> np.dtype:
    """
    Smallest SciPy-compatible index-pointer dtype for cumulative nnz counts.
    """
    return np.int32 if int(nnz) <= np.iinfo(np.int32).max else np.int64


def _rng_integers(
    rng: np.random.Generator,
    high: int,
    *,
    size,
    dtype,
) -> np.ndarray:
    """
    Generator.integers with dtype fallback for older NumPy versions.
    """
    try:
        return rng.integers(high, size=size, dtype=dtype)
    except TypeError:
        return rng.integers(high, size=size).astype(dtype, copy=False)


def _unique_rows_mask(values: np.ndarray) -> np.ndarray:
    """
    Return True for rows with no duplicate entries.
    """
    if values.shape[1] <= 1:
        return np.ones(values.shape[0], dtype=np.bool_)
    ordered = np.sort(values, axis=1)
    return np.all(ordered[:, 1:] != ordered[:, :-1], axis=1)


def _draw_unique_rows_by_rejection(
    rng: np.random.Generator,
    n_items: int,
    n_draw: int,
    n_rows: int,
    dtype,
) -> np.ndarray:
    """
    Draw n_rows ordered samples of length n_draw without replacement.

    IID draws conditioned on uniqueness are uniform over ordered samples without
    replacement. For sparse matrices where n_draw << n_items, rejection is
    typically much cheaper than one rng.choice call per entity.
    """
    out = _rng_integers(rng, n_items, size=(n_rows, n_draw), dtype=dtype)
    bad = ~_unique_rows_mask(out)
    while np.any(bad):
        n_bad = int(np.count_nonzero(bad))
        redraw = _rng_integers(rng, n_items, size=(n_bad, n_draw), dtype=dtype)
        out[bad] = redraw
        bad_idx = np.flatnonzero(bad)
        bad[bad_idx] = ~_unique_rows_mask(redraw)
    return out


def _fill_complement_indices(out: np.ndarray, n_items: int, omitted: np.ndarray) -> None:
    """
    Fill ``out`` with all item IDs except the sorted/unsorted omitted IDs.
    """
    omitted = np.sort(np.asarray(omitted, dtype=np.int64))
    pos = 0
    prev = 0
    dtype = out.dtype
    for miss_raw in omitted:
        miss = int(miss_raw)
        if miss > prev:
            n = miss - prev
            out[pos:pos + n] = np.arange(prev, miss, dtype=dtype)
            pos += n
        prev = miss + 1
    if prev < int(n_items):
        out[pos:] = np.arange(prev, int(n_items), dtype=dtype)


def _draw_omitted_rows(
    rng: np.random.Generator,
    n_items: int,
    n_omit: int,
    n_rows: int,
    dtype,
    *,
    small_degree_max: int,
) -> np.ndarray:
    """
    Draw omitted item IDs for complement sampling.
    """
    if n_omit <= 0:
        return np.empty((int(n_rows), 0), dtype=dtype)
    if n_omit == 1:
        return _rng_integers(rng, n_items, size=(int(n_rows), 1), dtype=dtype)
    if n_omit == 2:
        a = _rng_integers(rng, n_items, size=int(n_rows), dtype=dtype)
        b = _rng_integers(rng, n_items - 1, size=int(n_rows), dtype=dtype)
        b += (b >= a).astype(dtype, copy=False)
        return np.column_stack((a, b))

    collision_score = (float(n_omit) * float(n_omit - 1)) / (2.0 * float(n_items))
    if n_omit <= int(small_degree_max) or collision_score <= 0.25:
        return _draw_unique_rows_by_rejection(rng, n_items, n_omit, int(n_rows), dtype)

    out = np.empty((int(n_rows), int(n_omit)), dtype=dtype)
    for i in range(int(n_rows)):
        out[i] = rng.choice(n_items, size=int(n_omit), replace=False).astype(dtype, copy=False)
    return out


def _fill_fixed_degree_complement_blocks(
    entities: np.ndarray,
    n_items: int,
    n_omit: int,
    indptr: np.ndarray,
    indices: np.ndarray,
    rng: np.random.Generator,
    *,
    small_degree_max: int,
    target_temp_entries: int,
) -> None:
    """
    Fill equal-degree entities by sampling the omitted complement.
    """
    n_items = int(n_items)
    n_omit = int(n_omit)
    if entities.size == 0:
        return

    if n_omit == 0:
        all_items = np.arange(n_items, dtype=indices.dtype)
        for entity in entities:
            s = int(indptr[entity])
            e = int(indptr[entity + 1])
            indices[s:e] = all_items
        return

    block_size = max(1, int(target_temp_entries) // max(1, n_omit))
    for start in range(0, int(entities.size), block_size):
        ent = entities[start:start + block_size]
        omitted = _draw_omitted_rows(
            rng,
            n_items,
            n_omit,
            int(ent.size),
            indices.dtype,
            small_degree_max=small_degree_max,
        )
        for row_pos, entity in enumerate(ent):
            s = int(indptr[entity])
            e = int(indptr[entity + 1])
            _fill_complement_indices(indices[s:e], n_items, omitted[row_pos])


def _fill_random_key_subset_blocks(
    entities: np.ndarray,
    n_items: int,
    n_draw: int,
    indptr: np.ndarray,
    indices: np.ndarray,
    rng: np.random.Generator,
    *,
    sort_within_entity: bool,
    target_key_entries: int,
) -> None:
    """
    Fill equal-degree entities using random keys and argpartition.

    IID random keys induce a uniform random ordering of items per entity; taking
    the n_draw smallest keys gives a uniform subset without replacement.
    """
    n_items = int(n_items)
    n_draw = int(n_draw)
    if entities.size == 0:
        return

    dtype = indices.dtype
    if n_draw == n_items:
        all_items = np.arange(n_items, dtype=dtype)
        for entity in entities:
            s = int(indptr[entity])
            e = int(indptr[entity + 1])
            indices[s:e] = all_items
        return

    use_complement = n_draw > n_items // 2
    keep_n = n_items - n_draw if use_complement else n_draw
    block_size = max(1, int(target_key_entries) // n_items)
    arange_keep = np.arange(keep_n, dtype=indptr.dtype)

    keep_mask = None
    for start in range(0, int(entities.size), block_size):
        ent = entities[start:start + block_size]
        keys = rng.random((int(ent.size), n_items))
        chosen = np.argpartition(keys, kth=keep_n - 1, axis=1)[:, :keep_n]

        if use_complement:
            if keep_mask is None:
                keep_mask = np.empty(n_items, dtype=np.bool_)
            for row_pos, entity in enumerate(ent):
                s = int(indptr[entity])
                e = int(indptr[entity + 1])
                keep_mask.fill(True)
                keep_mask[chosen[row_pos]] = False
                indices[s:e] = np.flatnonzero(keep_mask).astype(dtype, copy=False)
            continue

        if sort_within_entity:
            chosen.sort(axis=1)

        offsets = indptr[ent, None] + arange_keep
        indices[offsets.ravel()] = chosen.ravel().astype(dtype, copy=False)


def _fill_fixed_degree_random_indices(
    degrees: np.ndarray,
    n_items: int,
    indptr: np.ndarray,
    indices: np.ndarray,
    rng: np.random.Generator,
    *,
    sort_within_entity: bool = False,
    small_degree_max: int = 16,
    target_temp_entries: int = 5_000_000,
    random_key_min_density: float = 0.01,
    target_key_entries: int = 20_000_000,
    debug_callback=None,
    debug_label: str = "fixed_degree",
) -> None:
    """
    Fill sparse indices for a fixed-degree null model.

    Entities with the same degree are sampled in blocks to collapse many tiny
    rng.choice calls into vectorised draws. Dense degree groups use random-key
    argpartition to avoid one rng.choice call per entity.
    """
    n_items = int(n_items)
    if n_items <= 0:
        return

    degrees = np.asarray(degrees, dtype=np.int64)
    nonzero = np.flatnonzero(degrees > 0)
    if nonzero.size == 0:
        return

    dtype = indices.dtype
    nz_degrees = degrees[nonzero]
    order = np.argsort(nz_degrees, kind="stable")
    grouped_entities = nonzero[order]
    grouped_degrees = nz_degrees[order]
    unique_degrees, starts = np.unique(grouped_degrees, return_index=True)
    ends = np.r_[starts[1:], grouped_entities.size]

    for k_raw, group_start, group_end in zip(unique_degrees, starts, ends):
        k = int(k_raw)
        entities = grouped_entities[int(group_start):int(group_end)]
        if entities.size == 0:
            continue
        if k > n_items:
            raise ValueError("Entity degree exceeds number of available items.")

        n_omit = n_items - k
        density = float(k) / float(n_items)
        complement_density = float(n_omit) / float(n_items)
        collision_score = (float(k) * float(k - 1)) / (2.0 * float(n_items)) if k > 1 else 0.0

        if k == n_items:
            _record_debug_event(
                debug_callback,
                sampler=debug_label,
                algorithm="full",
                degree=k,
                entities=int(entities.size),
                n_items=n_items,
                estimated_temp_entries=int(n_items),
            )
            _fill_fixed_degree_complement_blocks(
                entities,
                n_items,
                0,
                indptr,
                indices,
                rng,
                small_degree_max=small_degree_max,
                target_temp_entries=int(target_temp_entries),
            )
            continue

        use_complement = (
            k > n_items // 2
            and n_omit > 0
            and (n_omit <= int(small_degree_max) or complement_density <= 0.01)
        )
        if use_complement:
            _record_debug_event(
                debug_callback,
                sampler=debug_label,
                algorithm="complement_rejection",
                degree=k,
                omitted_degree=int(n_omit),
                entities=int(entities.size),
                n_items=n_items,
                estimated_temp_entries=min(int(target_temp_entries), int(entities.size) * int(n_omit)),
            )
            _fill_fixed_degree_complement_blocks(
                entities,
                n_items,
                n_omit,
                indptr,
                indices,
                rng,
                small_degree_max=small_degree_max,
                target_temp_entries=int(target_temp_entries),
            )
            continue

        use_rejection = k <= int(small_degree_max) or (density <= 0.01 and collision_score <= 0.25)
        use_random_key = (
            not use_rejection
            and min(k, n_omit) / n_items >= float(random_key_min_density)
            and n_items <= int(target_key_entries)
        )

        if use_random_key:
            _record_debug_event(
                debug_callback,
                sampler=debug_label,
                algorithm="random_key_argpartition",
                degree=k,
                entities=int(entities.size),
                n_items=n_items,
                estimated_temp_entries=min(
                    int(target_key_entries),
                    int(max(1, int(target_key_entries) // n_items)) * n_items,
                ),
            )
            _fill_random_key_subset_blocks(
                entities,
                n_items,
                k,
                indptr,
                indices,
                rng,
                sort_within_entity=sort_within_entity,
                target_key_entries=int(target_key_entries),
            )
            continue

        if use_rejection:
            algorithm = "rejection"
        elif n_omit < k:
            algorithm = "choice_complement_loop"
        else:
            algorithm = "choice_loop"
        _record_debug_event(
            debug_callback,
            sampler=debug_label,
            algorithm=algorithm,
            degree=k,
            entities=int(entities.size),
            n_items=n_items,
            estimated_temp_entries=min(int(target_temp_entries), int(entities.size) * max(1, min(k, n_omit))),
        )

        block_size = max(1, int(target_temp_entries) // max(1, k))
        for start in range(0, int(entities.size), block_size):
            ent = entities[start:start + block_size]
            m = int(ent.size)

            if k == 1:
                indices[indptr[ent]] = _rng_integers(rng, n_items, size=m, dtype=dtype)
                continue

            if k == 2:
                a = _rng_integers(rng, n_items, size=m, dtype=dtype)
                b = _rng_integers(rng, n_items - 1, size=m, dtype=dtype)
                b += (b >= a).astype(dtype, copy=False)
                if sort_within_entity:
                    lo = np.minimum(a, b)
                    hi = np.maximum(a, b)
                    vals = np.column_stack((lo, hi))
                else:
                    vals = np.column_stack((a, b))
            elif use_rejection:
                vals = _draw_unique_rows_by_rejection(rng, n_items, k, m, dtype)
                if sort_within_entity:
                    vals.sort(axis=1)
            elif n_omit < k:
                _fill_fixed_degree_complement_blocks(
                    ent,
                    n_items,
                    n_omit,
                    indptr,
                    indices,
                    rng,
                    small_degree_max=small_degree_max,
                    target_temp_entries=int(target_temp_entries),
                )
                continue
            else:
                for entity in ent:
                    s = int(indptr[entity])
                    e = int(indptr[entity + 1])
                    draw = rng.choice(n_items, size=k, replace=False)
                    if sort_within_entity:
                        draw.sort()
                    indices[s:e] = draw.astype(dtype, copy=False)
                continue

            offsets = indptr[ent, None] + np.arange(k, dtype=indptr.dtype)
            indices[offsets.ravel()] = vals.ravel()


def _check_linear_population_size(n_cells: int) -> None:
    """
    NumPy integer sampling for direct EE uses signed int64 cell IDs.
    """
    if int(n_cells) > np.iinfo(np.int64).max:
        raise ValueError("EE direct sampler requires n_rows * n_cols to fit in int64.")


def _sample_ee_linear_ids(
    rng: np.random.Generator,
    population_size: int,
    n_draw: int,
    *,
    debug_callback=None,
) -> np.ndarray:
    """
    Sample EE support cell IDs uniformly without replacement.
    """
    population_size = int(population_size)
    n_draw = int(n_draw)
    if n_draw <= 0:
        return np.empty(0, dtype=np.int64)
    if n_draw > population_size:
        raise ValueError("Cannot sample more unique IDs than the population size.")

    density = float(n_draw) / float(population_size) if population_size else 0.0
    _record_debug_event(
        debug_callback,
        sampler="EE",
        algorithm="numpy_choice",
        population_size=population_size,
        n_draw=n_draw,
        density=density,
        estimated_temp_entries=n_draw,
    )
    return rng.choice(population_size, size=n_draw, replace=False, shuffle=False).astype(np.int64, copy=False)


def _csr_from_linear_ids(
    linear_ids: np.ndarray,
    *,
    n_rows: int,
    n_cols: int,
) -> sp.csr_matrix:
    """
    Build a sorted CSR matrix from sampled linear cell IDs.
    """
    linear_ids = np.asarray(linear_ids, dtype=np.int64)
    linear_ids.sort()
    N = int(linear_ids.size)
    indptr = np.empty(int(n_rows) + 1, dtype=_sparse_indptr_dtype(N))
    indptr[0] = 0

    if N == 0:
        return sp.csr_matrix((int(n_rows), int(n_cols)), dtype=_OUT_DTYPE)

    rows = linear_ids // int(n_cols)
    row_counts = np.bincount(rows, minlength=int(n_rows)).astype(np.int64, copy=False)
    np.cumsum(row_counts, dtype=indptr.dtype, out=indptr[1:])
    cols = (linear_ids - rows * int(n_cols)).astype(_sparse_indices_dtype(n_cols), copy=False)
    data = np.ones(N, dtype=_OUT_DTYPE)
    return sp.csr_matrix((data, cols, indptr), shape=(int(n_rows), int(n_cols)))


def fe_fixed_rows_equiprob_cols(
    X_csr: sp.csr_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
    debug_callback=None,
) -> Iterable[sp.csr_matrix]:
    """
    FE model: fixed row totals, columns equiprobable without replacement per row.
    """
    n_rows, n_cols = X_csr.shape
    row_deg = (X_csr.indptr[1:] - X_csr.indptr[:-1]).astype(np.int64, copy=False)
    N = int(row_deg.sum())
    
    indptr = np.empty(n_rows + 1, dtype=_sparse_indptr_dtype(N))
    indptr[0] = 0
    np.cumsum(row_deg, dtype=indptr.dtype, out=indptr[1:])
    
    data = np.ones(N, dtype=_OUT_DTYPE)
    rng = _rng_from_state(random_state)
    
    for _ in range(int(n_reps)):
        indices = np.empty(N, dtype=_sparse_indices_dtype(n_cols))
        _fill_fixed_degree_random_indices(
            row_deg,
            n_cols,
            indptr,
            indices,
            rng,
            sort_within_entity=bool(sort_indices),
            debug_callback=debug_callback,
            debug_label="FE",
        )
        Y = sp.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
        if sort_indices:
            Y.sort_indices()
        yield Y


def ef_equiprob_rows_fixed_cols(
    X_csc: sp.csc_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
    debug_callback=None,
) -> Iterable[sp.csr_matrix]:
    """
    EF model: fixed column totals, rows equiprobable without replacement per column.
    Construct in CSC and convert to CSR per replicate.
    """
    n_rows, n_cols = X_csc.shape
    col_deg = (X_csc.indptr[1:] - X_csc.indptr[:-1]).astype(np.int64, copy=False)
    N = int(col_deg.sum())
    
    indptr = np.empty(n_cols + 1, dtype=_sparse_indptr_dtype(N))
    indptr[0] = 0
    np.cumsum(col_deg, dtype=indptr.dtype, out=indptr[1:])
    
    data = np.ones(N, dtype=_OUT_DTYPE)
    rng = _rng_from_state(random_state)
    
    for _ in range(int(n_reps)):
        indices = np.empty(N, dtype=_sparse_indices_dtype(n_rows))
        _fill_fixed_degree_random_indices(
            col_deg,
            n_rows,
            indptr,
            indices,
            rng,
            sort_within_entity=bool(sort_indices),
            debug_callback=debug_callback,
            debug_label="EF",
        )
        Y = sp.csc_matrix((data, indices, indptr), shape=(n_rows, n_cols)).tocsr()
        if sort_indices:
            Y.sort_indices()
        yield Y


def ee_equiprobable(
    X_csr: sp.csr_matrix,
    n_reps: int,
    random_state: Optional[int | np.random.Generator] = None,
    sort_indices: bool = False,
    debug_callback=None,
) -> Iterable[sp.csr_matrix]:
    """
    EE model: preserve total fill and sample occupied cells uniformly without replacement.
    """
    n_rows, n_cols = map(int, X_csr.shape)
    N = int(X_csr.nnz)
    n_cells = int(n_rows) * int(n_cols)
    _check_linear_population_size(n_cells)
    rng = _rng_from_state(random_state)
    
    if N == 0 or n_cells == 0:
        for _ in range(int(n_reps)):
            yield sp.csr_matrix((n_rows, n_cols), dtype=_OUT_DTYPE)
        return
        
    if N > n_cells:
        raise ValueError("nnz exceeds total number of cells; invalid input matrix.")

    for _ in range(int(n_reps)):
        ids = _sample_ee_linear_ids(
            rng,
            n_cells,
            N,
            debug_callback=debug_callback,
        )
        Y = _csr_from_linear_ids(ids, n_rows=n_rows, n_cols=n_cols)
        if sort_indices:
            Y.sort_indices()
        yield Y


# -----------------------------------------------------------------------------
# Samplers
# -----------------------------------------------------------------------------

class NullSampler(Protocol):
    """
    Sampler that yields CSR int8 presence matrices.
    """

    def sample(self, n_reps: int, *, seed: Optional[int] = None) -> Iterable[sp.csr_matrix]:
        """
        Yield n_reps null matrices in CSR format with int8 data.
        """
        ...


class DirectNullSampler:
    """
    Null sampler for FE/EF/EE models.
    Each call to sample() produces an independent replicate stream controlled by seed.
    """

    __slots__ = ("X", "model", "sort_indices", "default_random_state", "debug_callback", "ee_strategy")

    def __init__(
        self,
        X: sp.spmatrix,
        model: str,
        *,
        sort_indices: bool,
        default_random_state: Optional[int | np.random.Generator],
        debug_callback=None,
        ee_strategy: Optional[str] = None,
    ):
        self.X = X
        self.model = str(model).upper()
        self.sort_indices = bool(sort_indices)
        self.default_random_state = default_random_state
        self.debug_callback = debug_callback
        self.ee_strategy = ee_strategy

    def sample(self, n_reps: int, *, seed: Optional[int] = None) -> Iterable[sp.csr_matrix]:
        n_reps = int(n_reps)
        if n_reps <= 0:
            return iter(())

        eff_state = self.default_random_state if seed is None else int(seed)

        if self.model == "FE":
            return fe_fixed_rows_equiprob_cols(
                self.X,
                n_reps,
                random_state=eff_state,
                sort_indices=self.sort_indices,
                debug_callback=self.debug_callback,
            )

        if self.model == "EF":
            return ef_equiprob_rows_fixed_cols(
                self.X,
                n_reps,
                random_state=eff_state,
                sort_indices=self.sort_indices,
                debug_callback=self.debug_callback,
            )

        if self.model == "EE":
            return ee_equiprobable(
                self.X,
                n_reps,
                random_state=eff_state,
                sort_indices=self.sort_indices,
                debug_callback=self.debug_callback,
            )

        raise ValueError(f"Unknown or unsupported direct null model: {self.model}")


class FFCurveballSampler:
    """
    Stateful Curveball Markov chain sampler.
    Burn-in is executed once during initialisation.
    The first replicate is emitted immediately after burn-in.
    """

    __slots__ = (
        "X",
        "indptr",
        "indices",
        "nonempty",
        "ws",
        "rng",
        "burn_in_steps",
        "steps_per_rep",
        "sort_indices",
        "_first",
    )

    def __init__(
        self,
        X: sp.spmatrix,
        *,
        burn_in_steps: int,
        steps_per_rep: int,
        rng: np.random.Generator,
        sort_indices: bool,
        burn_q=None,
        burn_every: int = 0,
    ):
        self.X = X
        self.indptr = X.indptr
        self.indices = X.indices

        nnz_per_entity = self.indptr[1:] - self.indptr[:-1]
        self.nonempty = np.flatnonzero(nnz_per_entity > 0).astype(np.int64, copy=False)

        self.ws = _CurveballWorkspace(dtype=self.indices.dtype)
        self.rng = rng

        self.burn_in_steps = int(burn_in_steps)
        self.steps_per_rep = int(steps_per_rep)
        self.sort_indices = bool(sort_indices)
        self._first = True

        self._burn_in(burn_q=burn_q, burn_every=int(burn_every or 0))

    def _burn_in(self, *, burn_q=None, burn_every: int = 0) -> None:
        burn = int(self.burn_in_steps)
        if burn <= 0 or self.nonempty.size < 2:
            if burn_q is not None:
                burn_q.put(("DONE", 1))
            return

        if burn_q is None or burn_every <= 0:
            _curveball_steps(self.indptr, self.indices, self.nonempty, self.rng, burn, self.ws)
            if burn_q is not None:
                burn_q.put(("DONE", 1))
            return

        done = 0
        while done < burn:
            step = min(burn_every, burn - done)
            _curveball_steps(self.indptr, self.indices, self.nonempty, self.rng, step, self.ws)
            done += step
            burn_q.put(int(step))

        burn_q.put(("DONE", 1))

    def _snapshot(self) -> sp.csr_matrix:
        # Y = self.X.tocsr(copy=True)
        # Y.data = np.ones(Y.indices.size, dtype=_OUT_DTYPE)
        # if self.sort_indices:
            # Y.sort_indices()
        # return Y
        return self.X

    def sample(self, n_reps: int, *, seed: Optional[int] = None) -> Iterable[sp.csr_matrix]:
        n_reps = int(n_reps)
        if n_reps <= 0:
            return iter(())

        for _ in range(n_reps):
            if self._first:
                self._first = False
                yield self._snapshot()
                continue

            if self.steps_per_rep > 0:
                _curveball_steps(
                    self.indptr,
                    self.indices,
                    self.nonempty,
                    self.rng,
                    self.steps_per_rep,
                    self.ws,
                )
            yield self._snapshot()


def _ff_defaults(
    shape: tuple[int, int],
    burn_in_steps: Optional[int],
    steps_per_rep: Optional[int],
) -> tuple[int, int]:
    """
    Choose default FF chain parameters from matrix shape.
    """
    n_rows, n_cols = map(int, shape)
    basis = max(n_rows, n_cols)
    burn = int(burn_in_steps) if burn_in_steps is not None else max(1000, 5 * basis)
    steps = int(steps_per_rep) if steps_per_rep is not None else max(basis, 10)
    return burn, steps


def make_null_sampler(
    X: sp.spmatrix,
    model: SIMMODEL,
    *,
    random_state: Optional[int | np.random.Generator] = None,
    prepared: bool = False,
    copy: bool = True,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    sort_indices: bool = False,
    burn_q=None,
    burn_every: int = 0,
    debug_callback=None,
    ee_strategy: Optional[str] = None,
) -> NullSampler:
    """
    Construct a null sampler that yields CSR int8 presence matrices.

    FF is a stateful Curveball chain with one-time burn-in.
    FE/EF/EE are direct samplers producing independent replicates per seed.
    """
    model_u = str(model).upper()

    if model_u == "FF":
        n_rows, n_cols = map(int, X.shape)
        if not prepared:
            fmt = "csr" if n_rows <= n_cols else "csc"
            Xp = prepare_presence_matrix(X, fmt=fmt, copy=bool(copy))
        else:
            Xp = X.copy() if copy else X

        burn, steps = _ff_defaults(Xp.shape, burn_in_steps, steps_per_rep)
        rng = _rng_from_state(random_state)

        return FFCurveballSampler(
            Xp,
            burn_in_steps=burn,
            steps_per_rep=steps,
            rng=rng,
            sort_indices=sort_indices,
            burn_q=burn_q,
            burn_every=burn_every,
        )

    if model_u == "EF":
        Xp = prepare_presence_matrix(X, fmt="csc", copy=bool(copy)) if not prepared else (X.copy() if copy else X)
        return DirectNullSampler(
            Xp,
            model_u,
            sort_indices=sort_indices,
            default_random_state=random_state,
            debug_callback=debug_callback,
            ee_strategy=ee_strategy,
        )

    if model_u in ("FE", "EE"):
        Xp = prepare_presence_matrix(X, fmt="csr", copy=bool(copy)) if not prepared else (X.copy() if copy else X)
        return DirectNullSampler(
            Xp,
            model_u,
            sort_indices=sort_indices,
            default_random_state=random_state,
            debug_callback=debug_callback,
            ee_strategy=ee_strategy,
        )

    raise ValueError(f"Unknown or unsupported null model: {model_u}")


# -----------------------------------------------------------------------------
# Chunk planning and accumulator utilities
# -----------------------------------------------------------------------------

def _seed_seq_spawn(master_seed: Optional[int], n: int) -> list[int]:
    """
    Generate deterministic per-worker seeds using NumPy SeedSequence.
    """
    ss = np.random.SeedSequence(0 if master_seed is None else int(master_seed))
    return [int(s.generate_state(1, dtype=np.uint32)[0]) for s in ss.spawn(int(n))]


def _chunk_reps(n_reps: int, block: int) -> list[int]:
    """
    Split n_reps into chunks of size <= block.
    """
    out: list[int] = []
    n_reps = int(n_reps)
    block = max(1, int(block))
    while n_reps > 0:
        k = min(block, n_reps)
        out.append(int(k))
        n_reps -= k
    return out


def _merge_accumulators(a: dict, b: dict) -> dict:
    """
    Merge partial accumulator b into a in-place.
    """
    a["sum"] += b["sum"]
    a["sum2"] += b["sum2"]
    a["ge"] += b["ge"]
    a["n_ok"] += int(b["n_ok"])
    a["n_err"] += int(b["n_err"])
    a["n_done"] += int(b["n_done"])
    return a


def _finalise_accumulator(obs: np.ndarray, acc: dict) -> dict:
    """
    Convert summed moments into mean, sd, ses, and empirical p-values.
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
# Parallel reduction over null replicates
# -----------------------------------------------------------------------------

# Worker globals required by statistic functions and worker runner
_G_stat_fn = None
_G_obs = None

_G_term_cols = None
_G_nonterm_cols = None
_G_subset_idx = None
_G_n_rows = None
_G_N_T = None

_G_iA = None
_G_iB = None

_G_chunk_rows = None
_G_structure_do_nodf = None

_G_col_mask_T = None
_G_col_mask_notT = None
_G_sel_buf = None
_G_sel_buf_len = 0

_G_sampler = None

_G_burn_q = None
_G_burn_every = 0
_G_seed_q = None


def _row_counts_for_cols_csr(X_csr: sp.csr_matrix, col_mask: np.ndarray, out_len: int) -> np.ndarray:
    """
    Count, for each row, how many non-zeros fall within the selected columns.
    """
    indptr = X_csr.indptr
    sel = col_mask[X_csr.indices].astype(np.int64, copy=False)

    buf = _G_sel_buf
    buf[:sel.size] = sel
    buf[sel.size] = 0
    counts = np.add.reduceat(buf, indptr[:-1])

    if counts.size != out_len:
        counts = counts[:out_len]

    return counts


def _worker_init(
    X: sp.spmatrix,
    model: str,
    sort_indices: bool,
    burn_in_steps: Optional[int],
    steps_per_rep: Optional[int],
    stat_fn: Callable[[sp.csr_matrix], np.ndarray],
    obs: np.ndarray,
    **init_kwargs,
) -> None:
    """
    Cache shared inputs and preallocate buffers used by statistic functions.
    """
    global _G_stat_fn, _G_obs
    global _G_term_cols, _G_nonterm_cols, _G_subset_idx, _G_n_rows, _G_N_T
    global _G_iA, _G_iB
    global _G_chunk_rows, _G_structure_do_nodf
    global _G_col_mask_T, _G_col_mask_notT, _G_sel_buf, _G_sel_buf_len
    
    _G_stat_fn = stat_fn
    _G_obs = obs

    _G_term_cols = init_kwargs.get("term_cols", None)
    _G_nonterm_cols = init_kwargs.get("nonterm_cols", None)
    _G_subset_idx = init_kwargs.get("subset_idx", None)
    _G_n_rows = init_kwargs.get("n_rows", None)
    _G_N_T = init_kwargs.get("N_T", None)

    _G_iA = init_kwargs.get("iA", None)
    _G_iB = init_kwargs.get("iB", None)

    _G_chunk_rows = init_kwargs.get("chunk_rows", None)
    _G_structure_do_nodf = init_kwargs.get("structure_do_nodf", True)

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

    nnz = int(X.nnz)
    _G_sel_buf_len = nnz + 1
    _G_sel_buf = np.empty(_G_sel_buf_len, dtype=np.int64)


def _worker_init_wrap(
    X: sp.spmatrix,
    model: str,
    sort_indices: bool,
    burn_in_steps: Optional[int],
    steps_per_rep: Optional[int],
    stat_fn: Callable[[sp.csr_matrix], np.ndarray],
    obs: np.ndarray,
    init_kwargs: dict,
    burn_q,
    burn_every: int,
    seed_q,
) -> None:
    """
    Build the per-process null sampler and initialise statistic-function globals.
    """
    global _G_burn_q, _G_burn_every, _G_seed_q, _G_sampler
    
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
    
    model_u = str(model).upper()
    
    chain_seed = None
    if model_u == "FF" and _G_seed_q is not None:
        chain_seed = int(_G_seed_q.get())
        
    _G_sampler = make_null_sampler(
        X,
        model_u,
        random_state=(chain_seed if model_u == "FF" else None),
        prepared=True,
        copy=False,
        burn_in_steps=burn_in_steps,
        steps_per_rep=steps_per_rep,
        sort_indices=bool(sort_indices),
        burn_q=_G_burn_q if model_u == "FF" else None,
        burn_every=_G_burn_every if model_u == "FF" else 0,
    )


def _worker_run_chunk(task):
    """
    Evaluate statistic values over a chunk of null replicates and return partial sums.
    """
    n_reps_local, seed = task
    reps = int(n_reps_local)
    obs = _G_obs
    
    acc = {
        "sum": np.zeros_like(obs, dtype=float),
        "sum2": np.zeros_like(obs, dtype=float),
        "ge": np.zeros_like(obs, dtype=np.int64),
        "n_ok": 0,
        "n_err": 0,
        "n_done": 0,
    }
    
    for X_null in _G_sampler.sample(reps, seed=int(seed)):
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
            
    return acc, reps


# -----------------------------------------------------------------------------
# Statistic functions
# -----------------------------------------------------------------------------

def stat_fn_association_jaccard(X: sp.csr_matrix) -> np.ndarray:
    """
    Compute per-taxon association scores using counts in two column groups.
    """
    counts_T = _row_counts_for_cols_csr(X, _G_col_mask_T, _G_n_rows).astype(np.float64, copy=False)
    counts_notT = _row_counts_for_cols_csr(X, _G_col_mask_notT, _G_n_rows).astype(np.float64, copy=False)
    
    a_null = counts_T[_G_subset_idx]
    b_null = counts_notT[_G_subset_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(
            a_null,
            (b_null + _G_N_T),
            out=np.zeros_like(a_null, dtype=float),
            where=(b_null + _G_N_T) > 0,
        )


def stat_fn_cooccurrence_jaccard(X: sp.csr_matrix) -> np.ndarray:
    """
    Compute pairwise Jaccard overlap for a subset of taxa.
    """
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
    Compute [c_score, mean_jaccard, nodf] for a taxa × samples matrix.
    """
    from metacooc.structure import compute_c_score, mean_jaccard_dot, compute_nodf_streamed
    
    chunk_rows = int(_G_chunk_rows) if _G_chunk_rows is not None else 50_000
    do_nodf = bool(_G_structure_do_nodf)
    
    try:
        v0 = float(compute_c_score(X, chunk_rows=chunk_rows))
        v1 = float(mean_jaccard_dot(X, chunk_rows=chunk_rows))
        v2 = float(compute_nodf_streamed(X, chunk_rows=chunk_rows)) if do_nodf else np.nan
        
        if not (np.isfinite(v0) and np.isfinite(v1) and (np.isfinite(v2) or not do_nodf)):
            return None
            
        return np.array([v0, v1, v2], dtype=float)
        
    except MemoryError:
        return None
    except Exception:
        return None


# -----------------------------------------------------------------------------
# Parallel API
# -----------------------------------------------------------------------------

def parallel_null_reduce_vector(
    *,
    X: sp.spmatrix,
    model: str,
    n_reps: int,
    obs: np.ndarray,
    stat_fn: Callable[[sp.csr_matrix], np.ndarray],
    random_state: Optional[int] = None,
    n_workers: Optional[int] = None,
    chunksize: int = 1,
    sort_indices: bool = False,
    burn_in_steps: Optional[int] = None,
    steps_per_rep: Optional[int] = None,
    mp_start: str = "fork",
    progress_every: int = 1,
    **init_kwargs,
) -> dict:
    """
    Parallel reduction over null replicates for vector-valued statistics.
    
    FF uses an independent Markov chain per worker with one-time burn-in.
    FE/EF/EE use direct sampling with per-chunk deterministic seeds.
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
        
    # Default to a single worker unless explicitly requested
    if n_workers is None:
        n_workers = 1
    n_workers = max(1, int(n_workers))
   
   
    model_u = str(model).upper()
    
    chunk_n = int(progress_every) if progress_every and progress_every > 0 else 50
    chunk_n = max(1, chunk_n)
    
    chunk_sizes = _chunk_reps(n_reps, block=chunk_n)
    
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
    
    if model_u == "FF":
        n_rows, n_cols = map(int, X.shape)
        fmt = "csr" if n_rows <= n_cols else "csc"
        X = prepare_presence_matrix(X, fmt=fmt, copy=False)
        burn_in_steps, steps_per_rep = _ff_defaults((n_rows, n_cols), burn_in_steps, steps_per_rep)
        
        if n_rows <= n_cols:
            print("INFO: Curveball on rows (species); swapping samples between species")
        else:
            print("INFO: Curveball on columns (samples); swapping species between samples")
        print(f"INFO: Performing {int(burn_in_steps)} burn-in shuffles per chain before first statistics.")
        
    elif model_u == "EF":
        X = prepare_presence_matrix(X, fmt="csc", copy=False)
        
    elif model_u in ("FE", "EE"):
        X = prepare_presence_matrix(X, fmt="csr", copy=False)
        
    else:
        raise ValueError(f"Unknown or unsupported null model: {model_u}")
        
    import queue as _queue
    from tqdm import tqdm as _tqdm
    
    n_procs = min(n_workers, len(tasks))
    
    burn_q = None
    seed_q = None
    burn_every = 0
    
    if model_u == "FF":
        burn_q = ctx.Queue()
        seed_q = ctx.Queue()
        
        worker_seeds = _seed_seq_spawn(random_state, n_procs)
        for s in worker_seeds:
            seed_q.put(int(s))
            
        burn_every = 2000
        
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
            burn_q,
            burn_every,
            seed_q,
        ),
    ) as pool:
    
        if model_u == "FF":
            burn_steps = int(burn_in_steps or 0)
            if burn_steps > 0:
                total = burn_steps * n_procs
                done_workers = 0
                with _tqdm(total=total, desc=f"Burn-in (FF) - {n_procs} chains", dynamic_ncols=True) as pbar_burn:
                    while done_workers < n_procs:
                        try:
                            msg = burn_q.get(timeout=2)
                        except _queue.Empty:
                            continue

                        if isinstance(msg, tuple) and len(msg) == 2 and msg[0] == "DONE":
                            done_workers += int(msg[1])
                        else:
                            pbar_burn.update(int(msg))
                            
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
