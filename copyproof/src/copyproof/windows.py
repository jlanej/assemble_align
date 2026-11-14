"""Aggregate the positional k-mer track into fixed windows along every sequence.

For each window the tracked k-mers are split by their multiplicity in the whole release (m = 1,
2, ... up to a cap), and for each class we keep how many k-mer positions there are and the sum of
their read counts. Read counts divided by (positions x depth per copy) estimate how many copies of
that class of sequence the sample's genome holds; comparing it with m is the copy-number test.

Counts are winsorised per k-mer at the count expected for m + 4 copies, so a few k-mers from
collapsed high-copy repeats cannot dominate a window.

Neighbouring k-mers are not independent observations: two k-mers d bases apart are covered by
largely the same reads, so their counts covary by about depth x (1 - d / R), where R is the number
of k-mer start positions per read (read length - k + 1). For each window and class we therefore
also keep V, the sum of (1 - d / R) over all pairs of the class's k-mers in the window closer than
R (each k-mer paired with itself once). The variance of the class sum is then about
depth-per-k-mer x V, whatever the sampling rate, read length or clustering of the k-mers.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .track import Scan


def depth_per_copy(sc: Scan) -> dict:
    """Read depth per assembly copy (lambda), from the k-mers of the commonest multiplicity pattern.

    In a diploid assembly that pattern is "once in each haplotype"; its k-mers are overwhelmingly
    correct two-copy sequence, so the median of their read counts divided by 2 is the depth of one
    copy. This conditions on the assembly instead of fitting peaks of the read k-mer spectrum.
    """
    sp = sc.spectrum()
    mcols = [c for c in sp.columns if c.startswith("m_")]
    sp["m_tot"] = sp[mcols].sum(axis=1)
    by_pat = sp.groupby(mcols)["n_kmers"].sum().sort_values(ascending=False)
    pat = by_pat.index[0]
    pat = pat if isinstance(pat, tuple) else (pat,)
    sel = np.ones(len(sp), dtype=bool)
    for col, v in zip(mcols, pat):
        sel &= sp[col].values == v
    sub = sp[sel & (sp["count"] > 0)].sort_values("count")
    cum = sub["n_kmers"].cumsum().values
    med = float(sub["count"].values[np.searchsorted(cum, cum[-1] / 2)])
    m_tot = int(sum(pat))
    zero = float(sp.loc[sel & (sp["count"] == 0), "n_kmers"].sum())
    return {
        "lambda": med / m_tot,
        "pattern": dict(zip([c[2:] for c in mcols], map(int, pat))),
        "pattern_m": m_tot,
        "pattern_kmers": int(by_pat.iloc[0]),
        "pattern_median_count": med,
        "pattern_zero_fraction": zero / float(by_pat.iloc[0]),
    }


@dataclass
class Windows:
    """Per-window sufficient statistics for every sequence, concatenated."""

    W: int
    cap: int
    lam: float
    contigs: pd.DataFrame  # the scan's contig table plus win_start and n_win
    contig: np.ndarray  # contig index of each window
    start: np.ndarray  # window start on its contig
    n: np.ndarray  # (windows, cap + 1): tracked k-mer positions with total multiplicity m
    S: np.ndarray  # (windows, cap + 1): sum of their winsorised read counts
    n_sp: np.ndarray  # as n, for k-mers found only in this sequence's own group
    S_sp: np.ndarray
    zero: np.ndarray  # (windows, cap + 1): positions whose read count is 0
    V: np.ndarray  # (windows, cap + 1): variance multiplier of S (see the module docstring)
    n_all: np.ndarray  # tracked k-mer positions, including multiplicity above the cap
    gc: np.ndarray  # G+C fraction of the window's A/C/G/T bases
    acgt: np.ndarray  # A/C/G/T bases in the window
    # (windows, 3): sites (runs of adjacent k-mers) inherited from only the 1st parent, only the
    # 2nd parent, and found in neither parent
    hap: np.ndarray | None = None
    reff: int = 120  # k-mer start positions per read (read length - k + 1)

    def slice_of(self, ci: int) -> slice:
        row = self.contigs.iloc[ci]
        return slice(int(row.win_start), int(row.win_start) + int(row.n_win))

    def save(self, path: str) -> None:
        arrays = dict(
            W=self.W, cap=self.cap, lam=self.lam, contig=self.contig, start=self.start, n=self.n, S=self.S,
            n_sp=self.n_sp, S_sp=self.S_sp, zero=self.zero, V=self.V, n_all=self.n_all, gc=self.gc, acgt=self.acgt,
            reff=self.reff,
        )
        if self.hap is not None:
            arrays["hap"] = self.hap
        np.savez(path, **arrays)
        self.contigs.to_csv(str(path).removesuffix(".npz") + ".contigs.tsv", sep="\t", index=False)

    @staticmethod
    def load(path: str) -> "Windows":
        z = np.load(path)
        contigs = pd.read_csv(str(path).removesuffix(".npz") + ".contigs.tsv", sep="\t", dtype={"name": str, "group": str})
        return Windows(
            W=int(z["W"]), cap=int(z["cap"]), lam=float(z["lam"]), contigs=contigs, contig=z["contig"], start=z["start"],
            n=z["n"], S=z["S"], n_sp=z["n_sp"], S_sp=z["S_sp"], zero=z["zero"], V=z["V"], n_all=z["n_all"], gc=z["gc"],
            acgt=z["acgt"], hap=z["hap"] if "hap" in z.files else None, reff=int(z["reff"]),
        )


UNSUPPORTED_COLUMNS = ["contig", "start", "end", "unsupported", "tracked", "flank_copies", "flank_kmers", "single_copy"]


def count_cap(lam: float, cap: int) -> np.ndarray:
    """Winsorising ceiling of a k-mer's read count, by multiplicity: the count at m + 4 copies."""
    return np.ceil(lam * (np.arange(cap + 1) + 4.0)).astype(np.int64)


def read_length(sc: Scan) -> float:
    """Mean read length of the scan's first read set, from its database header."""
    db = sc.meta["reads"][0]["db"]
    return db["bases"] / max(db["records"], 1)


def pair_weights(key: np.ndarray, reff: int) -> np.ndarray:
    """For sorted integer positions `key`, the sum over neighbours within `reff` (itself included)
    of 1 - distance / reff."""
    n = len(key)
    if n == 0:
        return np.zeros(0)
    lo = np.searchsorted(key, key - reff, side="right")  # first neighbour with distance < reff
    hi = np.searchsorted(key, key + reff, side="left")  # one past the last such neighbour
    pre = np.concatenate([[0], np.cumsum(key, dtype=np.int64)])
    i = np.arange(n)
    right = (pre[hi] - pre[i + 1]) - key * (hi - i - 1)
    left = key * (i - lo) - (pre[i] - pre[lo])
    return (hi - lo) - (right + left) / float(reff)


def hapmer_masks(c: np.ndarray, child: int, parents: tuple[int, int], min_parent: int, min_child: int):
    """k-mers the child inherited from exactly one parent: present in that parent's reads
    (>= min_parent), absent from the other's (count 0), and present in the child's reads.
    The third mask marks k-mers present in the child's reads and in neither parent's."""
    a, b = c[:, parents[0]], c[:, parents[1]]
    ok = c[:, child] >= min_child
    return ok & (a >= min_parent) & (b == 0), ok & (b >= min_parent) & (a == 0), ok & (a == 0) & (b == 0)


def count_sites(pos: np.ndarray, win: np.ndarray, mask: np.ndarray, k: int, nw: int) -> np.ndarray:
    """Per window, the number of sites among the masked k-mers: k-mers closer than k to the
    previous one belong to the same site (one variant makes up to k adjacent k-mers)."""
    p = pos[mask]
    if len(p) == 0:
        return np.zeros(nw, dtype=np.int64)
    first = np.r_[True, np.diff(p) >= k]
    return np.bincount(win[mask][first], minlength=nw)


def unsupported_clusters(pos: np.ndarray, zero: np.ndarray, copies: np.ndarray, k: int, min_kmers: int, contig: str,
                         flank: int = 150, mult: np.ndarray | None = None) -> list[tuple]:
    """Runs of tracked k-mers with no read support, with the support of their surroundings.

    A wrong base or a false join in the assembly makes up to k unsupported k-mers in a row.
    Clusters are unsupported k-mers less than k apart. Each is returned as (contig, start, end,
    unsupported tracked k-mers, tracked k-mers in the span, median copies of the tracked k-mers
    within `flank` bases on either side, flank k-mers, unsupported k-mers that occur once in the
    release). `copies` is read count / (depth x multiplicity) per k-mer, so a flank value near 1
    means the cluster sits in normally covered sequence."""
    idx = np.flatnonzero(zero)
    if len(idx) == 0:
        return []
    p = pos[idx]
    brk = np.flatnonzero(np.diff(p) >= k) + 1
    starts = np.r_[0, brk]
    ends = np.r_[brk, len(p)]
    out = []
    for a, b in zip(starts, ends):
        if b - a < min_kmers:
            continue
        s0, s1 = int(p[a]), int(p[b - 1])
        i0, i1 = int(idx[a]), int(idx[b - 1]) + 1
        lo = int(np.searchsorted(pos, s0 - flank, side="left"))
        hi = int(np.searchsorted(pos, s1 + flank, side="right"))
        around = np.r_[copies[lo:i0], copies[i1:hi]]
        med = float(np.median(around)) if len(around) else float("nan")
        single = int((mult[idx[a:b]] == 1).sum()) if mult is not None else 0
        out.append((contig, s0, s1 + k, int(b - a), i1 - i0, round(med, 3), len(around), single))
    return out


def _contig_windows(prefix: str, ci: int, W: int, cap: int, lam: float, reff: int, parents, min_parent: int, min_child: int,
                    want_unsupported: bool, unsupported_min: int) -> dict:
    """Window statistics of one sequence (runs in a worker process)."""
    sc = _SCANS.get(prefix)
    if sc is None:
        sc = _SCANS[prefix] = Scan(prefix)
    ccap = count_cap(lam, cap)
    nb = W // sc.bin
    cols = cap + 1
    L = int(sc.contigs.at[ci, "length"])
    nw = -(-L // W)
    g = sc.group_index(ci)
    pos, m, c = sc.records(ci)
    win = pos // W
    mt = m.sum(axis=1, dtype=np.int64)
    elig = mt <= cap
    c0 = c[:, 0].astype(np.int64)
    idx = win[elig] * cols + mt[elig]
    cc = np.minimum(c0[elig], ccap[mt[elig]])
    size = nw * cols
    out = {"ci": ci, "nw": nw}
    out["n"] = np.bincount(idx, minlength=size).reshape(nw, cols).astype(np.int32)
    out["S"] = np.bincount(idx, weights=cc, minlength=size).reshape(nw, cols).astype(np.float32)
    sp = mt[elig] == m[elig, g]
    out["n_sp"] = np.bincount(idx[sp], minlength=size).reshape(nw, cols).astype(np.int32)
    out["S_sp"] = np.bincount(idx[sp], weights=cc[sp], minlength=size).reshape(nw, cols).astype(np.float32)
    out["zero"] = np.bincount(idx[c0[elig] == 0], minlength=size).reshape(nw, cols).astype(np.int32)
    out["unsupported"] = []
    if want_unsupported:
        cop = c0 / (lam * np.maximum(mt, 1))
        out["unsupported"] = unsupported_clusters(pos, (c0 == 0) & elig, cop, sc.k, unsupported_min, sc.contigs.at[ci, "name"],
                                                  mult=mt)
    # variance multiplier: pairs of same-class k-mers within a read length, inside one window
    V = np.zeros(size)
    pe, we, me = pos[elig], win[elig], mt[elig]
    key_all = pe + we * (W + reff)  # spreads windows apart so that pairs never straddle two
    for j in range(1, cols):
        sel = me == j
        if sel.any():
            V += np.bincount(idx[sel], weights=pair_weights(key_all[sel], reff), minlength=size)
    out["V"] = V.reshape(nw, cols).astype(np.float32)
    out["n_all"] = np.bincount(win, minlength=nw).astype(np.int32)
    gcb, acb = sc.composition(ci)
    pad = nw * nb - len(gcb)
    gcw = np.pad(gcb.astype(np.int64), (0, pad)).reshape(nw, nb).sum(axis=1)
    acw = np.pad(acb.astype(np.int64), (0, pad)).reshape(nw, nb).sum(axis=1)
    out["gc"] = (gcw / np.maximum(acw, 1)).astype(np.float32)
    out["acgt"] = acw.astype(np.int32)
    out["contig"] = np.full(nw, ci, dtype=np.int32)
    out["start"] = np.arange(nw, dtype=np.int64) * W
    if parents is not None:
        masks = hapmer_masks(c, 0, parents, min_parent, min_child)
        sites = [count_sites(pos, win, mk & elig, sc.k, nw) for mk in masks]
        out["hap"] = np.stack(sites, axis=1).astype(np.int32)
    return out


_SCANS: dict = {}


def compute_windows(
    sc: Scan,
    W: int = 1000,
    cap: int = 8,
    lam: float | None = None,
    parents: tuple[int, int] | None = None,
    min_parent: int = 3,
    min_child: int = 2,
    progress=None,
    unsupported: list | None = None,
    unsupported_min: int = 4,
    processes: int = 1,
) -> Windows:
    """Window statistics for the scan's first read set (the sample). With `processes` > 1 the
    sequences are handled in parallel (each worker holds one sequence's k-mers in memory)."""
    if not sc.reads:
        raise ValueError("the scan has no read set")
    if W % sc.bin:
        raise ValueError(f"window size must be a multiple of {sc.bin}")
    if lam is None:
        lam = depth_per_copy(sc)["lambda"]
    reff = max(int(round(read_length(sc))) - sc.k + 1, 1)
    args = (W, cap, float(lam), reff, parents, min_parent, min_child, unsupported is not None, unsupported_min)
    n_contigs = len(sc.contigs)
    results: list = [None] * n_contigs
    import multiprocessing as mp

    main_file = getattr(sys.modules.get("__main__"), "__file__", None)
    if processes > 1 and (main_file is None or not os.path.exists(main_file)):
        processes = 1  # worker processes re-import the main module, which a script read from stdin does not have
    if mp.parent_process() is not None:
        # Reached from inside a worker: the calling script lacks an `if __name__ == "__main__":`
        # guard, so every worker is re-running it. Stop instead of repeating the whole run.
        raise RuntimeError(
            "compute_windows was called from a worker process; guard the script's entry point with "
            "`if __name__ == \"__main__\":`"
        )
    if processes > 1 and n_contigs > 1:
        from concurrent.futures import ProcessPoolExecutor

        order = np.argsort(-sc.contigs["length"].values)  # longest first balances the workers
        with ProcessPoolExecutor(max_workers=processes, mp_context=mp.get_context("spawn")) as ex:
            futs = {ex.submit(_contig_windows, sc.prefix, int(ci), *args): int(ci) for ci in order}
            done = 0
            for fut, ci in futs.items():
                results[ci] = fut.result()
                done += 1
                if progress:
                    progress(done - 1, n_contigs)
    else:
        _SCANS[sc.prefix] = sc
        for ci in range(n_contigs):
            results[ci] = _contig_windows(sc.prefix, ci, *args)
            if progress:
                progress(ci, n_contigs)
    keys = ["contig", "start", "n", "S", "n_sp", "S_sp", "zero", "V", "n_all", "gc", "acgt"] + (["hap"] if parents is not None else [])
    out = {k: [r[k] for r in results] for k in keys}
    n_win = [r["nw"] for r in results]
    win_start = np.concatenate([[0], np.cumsum(n_win)[:-1]]).astype(int).tolist()
    if unsupported is not None:
        for r in results:
            unsupported.extend(r["unsupported"])
    contigs = sc.contigs.copy()
    contigs["win_start"] = win_start
    contigs["n_win"] = n_win
    cat = {k: np.concatenate(v) for k, v in out.items() if v}
    return Windows(
        W=W, cap=cap, lam=float(lam), contigs=contigs, contig=cat["contig"], start=cat["start"], n=cat["n"], S=cat["S"],
        n_sp=cat["n_sp"], S_sp=cat["S_sp"], zero=cat["zero"], V=cat["V"], n_all=cat["n_all"], gc=cat["gc"],
        acgt=cat["acgt"], hap=cat.get("hap"), reff=reff,
    )
