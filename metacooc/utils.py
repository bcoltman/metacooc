#!/usr/bin/env python3

import numpy as np
import scipy.sparse as sp
from typing import Iterable, Tuple

_RANK_PREFIXES = {
    "domain": "d__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}

_PREFIX_TO_RANK = {v: k for k, v in _RANK_PREFIXES.items()}

_RANK_ORDER = ["root", "domain", "phylum", "class", "order", "family", "genus", "species"]

def _parse_tokens(s: str) -> list[str]:
    # No regex, single strip per token
    parts = s.split(';')
    out: list[str] = []
    for t in parts:
        t = t.strip()
        if t:
            out.append(t)
    return out

def _token_rank(token: str) -> str | None:
    # assume already stripped
    if token.lower() == "root":
        return "root"
    if len(token) >= 3 and token[1:3] == "__":
        pref = token[:3]
        return _PREFIX_TO_RANK.get(pref)
    return None

def _terminal_rank_prefix(tokens: list[str]) -> str | None:
    if not tokens:
        return None
    last = tokens[-1]
    if len(last) >= 3 and last[1:3] == "__":
        pref = last[:3]
        if pref in _RANK_PREFIXES.values():  # or precompute a set of valid prefixes
            return pref
    return None

def _deepest_rank_token(search_string: str) -> tuple[str | None, str | None]:
    # Walk from right to left and return the deepest ranked token (ignores 'Root')
    for tok in reversed(_parse_tokens(search_string)):
        r = _token_rank(tok)
        if r is not None:
            return r, tok
    return None, None

def stream_csr_upper_threshold(
    M: sp.csr_matrix,
    threshold: float = 0.0,
    chunk_rows: int = 50_000,
) -> Iterable[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Stream strict upper-triangle entries (i < j) from CSR matrix M with val > threshold.
    Avoids COO conversion (no coo.row allocation).
    
    Yields arrays: (i, j, val)
    """
    if not sp.isspmatrix_csr(M):
        M = M.tocsr(copy=False)
    
    if M.shape[0] != M.shape[1]:
        raise ValueError("stream_csr_upper_threshold expects a square matrix")
        
    indptr = M.indptr
    indices = M.indices
    data = M.data
    n = M.shape[0]
    
    rows_out, cols_out, vals_out = [], [], []
    
    for r0 in range(0, n, chunk_rows):
        r1 = min(r0 + chunk_rows, n)
        
        rows_out.clear(); cols_out.clear(); vals_out.clear()
        
        for i in range(r0, r1):
            s, e = indptr[i], indptr[i + 1]
            if s == e:
                continue
                
            cols = indices[s:e]
            vals = data[s:e]
            
            m = (cols > i) & (vals > threshold)
            if np.any(m):
                c = cols[m]
                v = vals[m]
                rows_out.append(np.full(c.size, i, dtype=cols.dtype))
                cols_out.append(c)
                vals_out.append(v)
                
        if rows_out:
            yield (np.concatenate(rows_out),
                   np.concatenate(cols_out),
                   np.concatenate(vals_out))

def stream_edges(
    M_csr: sp.csr_matrix,
    threshold: float,
) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """
    Stream unique edges (iA < iB) from matrix.
    """
    for rows, cols, _ in stream_csr_upper_threshold(M_csr, threshold=threshold):
        # Here rows are iA, cols are iB already satisfying cols > rows
        yield rows, cols