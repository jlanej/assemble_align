"""Ground truth from a truth genome's k-mers.

When a complete, correct assembly of the same individual exists (T2T-HG002 Q100), the k-mers of
that truth genome say how many copies of each k-mer the genome really holds. Scanning a test
assembly with the truth's k-mer multiplicities as an extra "read set" then labels every tracked
k-mer of the test assembly: truth copies minus assembly copies. No alignment is involved, so the
labels are independent of any aligner's choices in repeats.

    cpk count --exact true --out truth  truth.fa
    cpk scan --asm ... --reads sample=reads.cpkdb --reads truth=truth.cpkdb ...
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .track import Scan


@dataclass
class TruthWindows:
    W: int
    cap: int
    contig: np.ndarray
    start: np.ndarray
    n: np.ndarray  # eligible tracked k-mer positions (assembly multiplicity <= cap)
    gain: np.ndarray  # positions where the truth holds more copies than the assembly
    loss: np.ndarray  # positions where the truth holds fewer copies (0 included)
    absent: np.ndarray  # positions whose k-mer is not in the truth at all
    net: np.ndarray  # sum over positions of (truth copies - assembly copies), each clipped to +-3


def truth_windows(sc: Scan, truth_index: int, W: int = 1000, cap: int = 8) -> TruthWindows:
    """Per-window truth labels along the scanned assembly."""
    out = {k: [] for k in ("contig", "start", "n", "gain", "loss", "absent", "net")}
    for ci in range(len(sc.contigs)):
        L = int(sc.contigs.at[ci, "length"])
        nw = -(-L // W)
        pos, m, c = sc.records(ci)
        mt = m.sum(axis=1, dtype=np.int64)
        elig = mt <= cap
        win = pos[elig] // W
        t = c[elig, truth_index].astype(np.int64)
        d = t - mt[elig]
        out["contig"].append(np.full(nw, ci, dtype=np.int32))
        out["start"].append(np.arange(nw, dtype=np.int64) * W)
        out["n"].append(np.bincount(win, minlength=nw).astype(np.int32))
        out["gain"].append(np.bincount(win[d > 0], minlength=nw).astype(np.int32))
        out["loss"].append(np.bincount(win[d < 0], minlength=nw).astype(np.int32))
        out["absent"].append(np.bincount(win[t == 0], minlength=nw).astype(np.int32))
        out["net"].append(np.bincount(win, weights=np.clip(d, -3, 3), minlength=nw).astype(np.float32))
    cat = {k: np.concatenate(v) for k, v in out.items()}
    return TruthWindows(W=W, cap=cap, **cat)


def truth_segments(tw: TruthWindows, contigs: pd.DataFrame, min_frac: float = 0.5, min_kmers: int = 30, min_win: int = 2,
                   merge_gap: int = 1) -> pd.DataFrame:
    """Runs of windows where most eligible k-mers have more (gain) or fewer (loss) copies in the
    truth than in the assembly."""
    rows = []
    offset = 0
    for ci in range(len(contigs)):
        nw = int((tw.contig == ci).sum())
        sl = slice(offset, offset + nw)
        offset += nw
        n = tw.n[sl].astype(float)
        name = contigs.at[ci, "name"]
        L = int(contigs.at[ci, "length"])
        for direction, arr in (("gain", tw.gain[sl]), ("loss", tw.loss[sl])):
            flag = (n >= min_kmers) & (arr >= min_frac * n)
            idx = np.flatnonzero(flag)
            if len(idx) == 0:
                continue
            brk = np.flatnonzero(np.diff(idx) > merge_gap + 1) + 1
            for run in np.split(idx, brk):
                a, b = int(run[0]), int(run[-1]) + 1
                if len(run) < min_win:
                    continue
                kk = n[a:b].sum()
                rows.append(
                    {
                        "contig": name, "start": a * tw.W, "end": min(b * tw.W, L), "direction": direction,
                        "n_win": len(run), "kmers": int(kk), "frac": round(float(arr[a:b].sum() / max(kk, 1)), 3),
                        "mean_delta": round(float(tw.net[sl][a:b].sum() / max(kk, 1)), 3),
                        "absent_frac": round(float(tw.absent[sl][a:b].sum() / max(kk, 1)), 3),
                    }
                )
    return pd.DataFrame(rows, columns=["contig", "start", "end", "direction", "n_win", "kmers", "frac", "mean_delta", "absent_frac"])


def match_segments(calls: pd.DataFrame, truth: pd.DataFrame, min_overlap: float = 0.5):
    """Mark calls supported by a truth segment of the same direction, and truth segments found by
    a call. Overlap is measured as a fraction of the call (precision side) or of the truth segment
    (recall side)."""
    def by_contig(df):
        return {c: g.sort_values("start") for c, g in df.groupby("contig")}

    tb, cb = by_contig(truth), by_contig(calls)

    def covered(row, others):
        g = others.get(row.contig)
        if g is None:
            return 0.0
        g = g[(g["direction"] == row.direction) & (g["start"] < row.end) & (g["end"] > row.start)]
        if not len(g):
            return 0.0
        ov = (np.minimum(g["end"].values, row.end) - np.maximum(g["start"].values, row.start)).sum()
        return float(ov) / max(row.end - row.start, 1)

    calls = calls.copy()
    truth = truth.copy()
    calls["truth_overlap"] = [covered(r, tb) for r in calls.itertuples(index=False)]
    truth["call_overlap"] = [covered(r, cb) for r in truth.itertuples(index=False)]
    calls["true"] = calls["truth_overlap"] >= min_overlap
    truth["found"] = truth["call_overlap"] >= min_overlap
    return calls, truth
