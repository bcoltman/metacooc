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
    chunk_rows: int = 10_000,
) -> Iterable[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Stream strict upper-triangle entries (i < j) from a square CSR matrix
    with val > threshold, without converting to COO.
    
    Parameters
    ----------
    M : sp.csr_matrix
        Input square sparse matrix.
    threshold : float, default 0.0
        Only yield entries with value > threshold.
    chunk_rows : int, default 10_000
        Number of rows to process per block.
        
    Yields
    ------
    rows, cols, vals : np.ndarray
        Parallel arrays for one block of retained entries.
    """
    if not sp.isspmatrix_csr(M):
        M = M.tocsr(copy=False)
        
    if M.shape[0] != M.shape[1]:
        raise ValueError("stream_csr_upper_threshold expects a square matrix")
        
    if chunk_rows <= 0:
        raise ValueError("chunk_rows must be a positive integer")
        
    indptr = M.indptr
    indices = M.indices
    data = M.data
    n = M.shape[0]
    
    use_value_threshold = threshold > 0
    
    for r0 in range(0, n, chunk_rows):
        r1 = min(r0 + chunk_rows, n)
        
        total_keep = 0
        
        # pass 1: count surviving entries in this block
        for i in range(r0, r1):
            s = indptr[i]
            e = indptr[i + 1]
            if s == e:
                continue
                
            cols = indices[s:e]
            
            if use_value_threshold:
                vals = data[s:e]
                keep = (cols > i) & (vals > threshold)
            else:
                keep = (cols > i)
                
            total_keep += int(np.count_nonzero(keep))
            
        if total_keep == 0:
            continue
            
        rows_out = np.empty(total_keep, dtype=indices.dtype)
        cols_out = np.empty(total_keep, dtype=indices.dtype)
        vals_out = np.empty(total_keep, dtype=data.dtype)
        
        # pass 2: fill arrays
        pos = 0
        for i in range(r0, r1):
            s = indptr[i]
            e = indptr[i + 1]
            if s == e:
                continue
                
            cols = indices[s:e]
            vals = data[s:e]
            
            if use_value_threshold:
                keep = (cols > i) & (vals > threshold)
            else:
                keep = (cols > i)
                
            k = int(np.count_nonzero(keep))
            if k == 0:
                continue
                
            rows_out[pos:pos + k] = i
            cols_out[pos:pos + k] = cols[keep]
            vals_out[pos:pos + k] = vals[keep]
            pos += k
            
        yield rows_out, cols_out, vals_out


def stream_edge_values(
    M_csr: sp.csr_matrix,
    min_value: int,
    chunk_rows: int = 10_000,
) -> Iterable[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Stream strict upper-triangle entries (i < j) and their values from a sparse matrix.
    
    Returns
    -------
    rows, cols, vals
        Parallel arrays where vals[k] = M_csr[rows[k], cols[k]].
    """
    for rows, cols, vals in stream_csr_upper_threshold(
        M_csr,
        threshold=min_value,
        chunk_rows=chunk_rows,
    ):
        yield rows, cols, vals

def _stream_csr_entries(
    M: sp.csr_matrix,
    min_value: int = 0,
    chunk_rows: int = 10_000,
):
    """
    Stream all non-zero entries from any CSR matrix, including rectangular ones.
    Yields row indices, column indices, and values.
    """
    if not sp.isspmatrix_csr(M):
        M = M.tocsr(copy=False)
        
    indptr = M.indptr
    indices = M.indices
    data = M.data
    n_rows = M.shape[0]
    use_value_threshold = min_value > 0
    
    for r0 in range(0, n_rows, chunk_rows):
        r1 = min(r0 + chunk_rows, n_rows)
        
        rows_parts = []
        cols_parts = []
        vals_parts = []
        
        for i in range(r0, r1):
            s, e = indptr[i], indptr[i + 1]
            if s == e:
                continue
                
            cols = indices[s:e]
            vals = data[s:e]
            
            if use_value_threshold:
                keep = vals > min_value
                if not np.any(keep):
                    continue
                cols = cols[keep]
                vals = vals[keep]
                
            rows_parts.append(np.full(cols.size, i, dtype=np.int64))
            cols_parts.append(cols.astype(np.int64, copy=False))
            vals_parts.append(vals)
            
        if rows_parts:
            yield (
                np.concatenate(rows_parts),
                np.concatenate(cols_parts),
                np.concatenate(vals_parts),
            )