"""Interval sets on named sequences (BED-like, 0-based half-open), built on sorted numpy arrays."""

from __future__ import annotations

import gzip
from collections import defaultdict

import numpy as np


class IntervalSet:
    """Merged, sorted intervals per sequence."""

    def __init__(self, by_contig: dict[str, np.ndarray] | None = None):
        self.by = {}
        for c, arr in (by_contig or {}).items():
            a = np.asarray(arr, dtype=np.int64).reshape(-1, 2)
            self.by[c] = _merge(a)

    @staticmethod
    def from_bed(path: str, min_cols: int = 3, keep=None) -> "IntervalSet":
        """Read a BED file; `keep(fields) -> bool` optionally filters lines."""
        acc = defaultdict(list)
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f:
            for line in f:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                t = line.rstrip("\n").split("\t")
                if len(t) < min_cols:
                    continue
                if keep is not None and not keep(t):
                    continue
                acc[t[0]].append((int(t[1]), int(t[2])))
        return IntervalSet({c: np.array(v, dtype=np.int64) for c, v in acc.items()})

    @staticmethod
    def from_rows(rows) -> "IntervalSet":
        acc = defaultdict(list)
        for c, s, e in rows:
            acc[c].append((int(s), int(e)))
        return IntervalSet({c: np.array(v, dtype=np.int64) for c, v in acc.items()})

    def total(self) -> int:
        return int(sum((a[:, 1] - a[:, 0]).sum() for a in self.by.values()))

    def __len__(self) -> int:
        return int(sum(len(a) for a in self.by.values()))

    def overlap(self, contig: str, start: int, end: int) -> int:
        """Bases of [start, end) covered by the set."""
        a = self.by.get(contig)
        if a is None or end <= start:
            return 0
        i = np.searchsorted(a[:, 1], start, side="right")
        j = np.searchsorted(a[:, 0], end, side="left")
        if j <= i:
            return 0
        s = np.maximum(a[i:j, 0], start)
        e = np.minimum(a[i:j, 1], end)
        return int((e - s).sum())

    def hits(self, contig: str, start: int, end: int) -> np.ndarray:
        """Intervals of the set that overlap [start, end)."""
        a = self.by.get(contig)
        if a is None:
            return np.zeros((0, 2), dtype=np.int64)
        i = np.searchsorted(a[:, 1], start, side="right")
        j = np.searchsorted(a[:, 0], end, side="left")
        return a[i:j]

    def union(self, other: "IntervalSet") -> "IntervalSet":
        out = {}
        for c in set(self.by) | set(other.by):
            parts = [x for x in (self.by.get(c), other.by.get(c)) if x is not None]
            out[c] = np.concatenate(parts)
        return IntervalSet(out)

    def pad(self, n: int) -> "IntervalSet":
        return IntervalSet({c: np.stack([np.maximum(a[:, 0] - n, 0), a[:, 1] + n], axis=1) for c, a in self.by.items()})

    def window_fraction(self, contig: str, length: int, W: int) -> np.ndarray:
        """Fraction of each W-sized window of a sequence covered by the set."""
        nw = -(-length // W)
        cov = np.zeros(nw + 1)
        a = self.by.get(contig)
        if a is None:
            return cov[:nw]
        # difference array over bases would be large; accumulate per window instead
        for s, e in a:
            s, e = max(int(s), 0), min(int(e), length)
            if e <= s:
                continue
            ws, we = s // W, (e - 1) // W
            if ws == we:
                cov[ws] += e - s
            else:
                cov[ws] += (ws + 1) * W - s
                cov[we] += e - we * W
                if we > ws + 1:
                    cov[ws + 1 : we] += W
        sizes = np.full(nw, W, dtype=float)
        sizes[-1] = length - (nw - 1) * W
        return cov[:nw] / sizes

    def to_rows(self):
        for c in sorted(self.by):
            for s, e in self.by[c]:
                yield c, int(s), int(e)


def _merge(a: np.ndarray) -> np.ndarray:
    if len(a) == 0:
        return a.reshape(0, 2)
    a = a[np.argsort(a[:, 0], kind="stable")]
    end = np.maximum.accumulate(a[:, 1])
    new = np.r_[True, a[1:, 0] > end[:-1]]
    idx = np.flatnonzero(new)
    starts = a[idx, 0]
    ends = np.r_[end[idx[1:] - 1], end[-1]]
    return np.stack([starts, ends], axis=1)
