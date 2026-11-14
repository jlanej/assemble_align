"""The other copies of a copy-number segment, and what they say about it.

A segment whose k-mers all follow one pattern is ambiguous: where the release holds a sequence
three times and the reads show two copies, "one copy too many" and "two thirds of the depth" are
the same observation. An extra copy that is identical to another copy (the usual false
duplication) has no k-mers of its own to decide it.

But its k-mers occur elsewhere in the release, and the two explanations differ in what they say
about those other copies. Whatever lowers sequencing depth on a sequence (base composition, a
structure the chemistry reads poorly) lowers it on every copy of that sequence, so the k-mers
unique to another copy, lying between the shared ones, must be depleted by the same factor.
Under a surplus copy they are read at full depth. So:

1. `cpk locate --regions` lists where the tracked k-mers of each segment occur in the release;
2. hits on one diagonal form a copy (`copy_loci`);
3. the one-copy k-mers of the copies' windows are tested at one copy against the segment's
   depth ratio (`resolve`), and their log-likelihood is added to the segment's own.

The copies are also written to the segment table, so an issue can name the sequences involved.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .cn import CnParams
from .engine import engine_path
from .model import Model
from .windows import Windows


def locate_segments(segs: pd.DataFrame, asm_args: list, k: int, sample: int, cap: int, prefix: Path, threads: int, log=None) -> pd.DataFrame:
    """Where the tracked k-mers of each segment occur in the release (the segment itself included).
    The query name is the segment's index."""
    bed = Path(f"{prefix}.regions.bed")
    hits = Path(f"{prefix}.hits.tsv")
    text = "".join(f"{r.contig}\t{int(r.start)}\t{int(r.end)}\t{i}\n" for i, r in zip(segs.index, segs.itertuples(index=False)))
    if not (bed.is_file() and hits.is_file() and bed.read_text() == text):
        hits.unlink(missing_ok=True)
        bed.write_text(text)
        subprocess.run([engine_path(), "locate", "--regions", str(bed), "--out", str(hits), "--k", str(k), "--sample", str(sample),
                        "--max-occ", str(cap), "--threads", str(threads)] + [str(a) for a in asm_args], check=True, stderr=log)
    return pd.read_csv(hits, sep="\t", dtype={"query": np.int64, "contig": str, "strand": str})


def copy_loci(h: pd.DataFrame, contig: str, start: int, end: int, k: int, min_kmers: int = 20, min_frac: float = 0.2,
              diag_tol: int = 2000) -> list[dict]:
    """Copies of one segment elsewhere in the release, from the located k-mers `h` of that segment.

    Hits of one sequence and strand that lie on one diagonal (within `diag_tol`, which allows for
    insertions and deletions between copies) are one copy. A copy must hold `min_kmers` of the
    segment's k-mers and `min_frac` of those that occur anywhere else."""
    pos = h["pos"].values
    self_hit = (h["contig"].values == contig) & (pos >= start) & (pos < end)
    h = h[~self_hit]
    if not len(h):
        return []
    n_multi = h["qpos"].nunique()
    need = max(min_kmers, min_frac * n_multi)
    out = []
    for (c, strand), g in h.groupby(["contig", "strand"], sort=False):
        if len(g) < need:
            continue
        q = g["qpos"].values.astype(np.int64)
        x = g["pos"].values.astype(np.int64)
        diag = x - q if strand == "+" else x + q
        order = np.argsort(diag, kind="stable")
        cut = np.flatnonzero(np.diff(diag[order]) > diag_tol) + 1
        for idx in np.split(order, cut):
            n = len(np.unique(q[idx]))
            if n < need:
                continue
            out.append({
                "contig": c, "start": int(x[idx].min()), "end": int(x[idx].max()) + k, "strand": strand, "kmers": int(n),
                "share": round(n / n_multi, 3), "q_start": int(q[idx].min()), "q_end": int(q[idx].max()) + k,
            })
    out.sort(key=lambda d: -d["kmers"])
    return out


def _windows_inside(w: Windows, index: dict, loci: list[dict]) -> np.ndarray:
    """Global indices of the windows that lie wholly inside any of the loci."""
    keep = []
    for d in loci:
        ci = index.get(d["contig"])
        if ci is None:
            continue
        sl = w.slice_of(ci)
        a = -(-d["start"] // w.W)
        b = d["end"] // w.W
        if b > a:
            keep.append(np.arange(sl.start + a, min(sl.start + b, sl.stop)))
    return np.unique(np.concatenate(keep)) if keep else np.zeros(0, dtype=np.int64)


def unique_at_copies(w: Windows, model: Model, p: CnParams, ok: np.ndarray, wins: np.ndarray, ratio: float, rel_ratio: float) -> dict:
    """The one-copy k-mers of the copies' windows, tested at one copy and at `ratio` copies.

    `ratio` is the depth of the segment's commonest class relative to its expectation, known to
    a relative `rel_ratio`. Returns the k-mers, their copies, both z-scores and the log-likelihood
    of "one copy" over "depleted like the segment" (positive: the copies are read at full depth)."""
    wins = wins[ok[wins]] if len(wins) else wins
    res = {"n": 0, "copies": float("nan"), "z_one": float("nan"), "z_ratio": float("nan"), "llr": 0.0, "windows": int(len(wins))}
    if not len(wins):
        return res
    n = w.n[wins, 1].astype(float)
    f = model.gc_factor(w.gc[wins].astype(float))
    N, F, FV, S = n.sum(), (f * n).sum(), (f * w.V[wins, 1].astype(float)).sum(), float(w.S[wins, 1].sum())
    res["n"] = int(N)
    if N < 1:
        return res
    om = model.omega_vec(w.cap)[1]
    lam = model.lam

    def fit(copies, rel):
        mu = lam * (copies + model.eps)
        E = mu * F
        sd = np.sqrt(om * mu * FV + (rel * E) ** 2 + 1.0)
        z = (S - E) / sd
        return -np.log(sd) - 0.5 * (p.nu + 1.0) * np.log1p(z * z / p.nu), z

    l1, z1 = fit(1.0, model.psi)
    lr, zr = fit(ratio, float(np.hypot(p.psi_local, rel_ratio)))
    res.update(copies=round(S / (lam * F), 3), z_one=round(float(z1), 2), z_ratio=round(float(zr), 2), llr=round(float(l1 - lr), 2))
    return res


def _ratio_of(w: Windows, model: Model, r, ok: np.ndarray) -> tuple[float, float]:
    """Depth of the segment's commonest class relative to its expectation, and the relative
    sampling error of that ratio."""
    m = int(r.dominant_m)
    wins = np.arange(int(r.w0), int(r.w1))
    wins = wins[ok[wins]]
    n = w.n[wins, m].astype(float)
    f = model.gc_factor(w.gc[wins].astype(float))
    F, FV, S = (f * n).sum(), (f * w.V[wins, m].astype(float)).sum(), float(w.S[wins, m].sum())
    if F <= 0 or S <= 0:
        return float("nan"), float("nan")
    om = model.omega_vec(w.cap)[m]
    return S / (model.lam * m * F), float(np.sqrt(om * (S / F) * FV + 1.0) / S)


def resolve(segs: pd.DataFrame, w: Windows, model: Model, p: CnParams, ok: np.ndarray, hits: pd.DataFrame, k: int,
            min_unique: int | None = None) -> pd.DataFrame:
    """Add the copies of each located segment to the table, and settle the ambiguous ones.

    New columns: `copies` (the other copies, as sequence:start-end:shared k-mers), `copy_unique`
    (one-copy k-mers at those copies and the copies the reads give them), `copy_llr` (their
    log-likelihood for full depth over the segment's depth ratio) and `resolved` ("copies" where
    this evidence changed the segment's kind)."""
    segs = segs.copy()
    for col, default in (("copies", ""), ("copy_unique", ""), ("copy_llr", 0.0), ("resolved", "")):
        segs[col] = default
    if not len(hits):
        return segs
    min_unique = p.min_seg_class if min_unique is None else min_unique
    index = {name: i for i, name in enumerate(w.contigs["name"])}
    by = {q: g for q, g in hits.groupby("query", sort=False)}
    for i, r in zip(segs.index, segs.itertuples(index=False)):
        h = by.get(i)
        if h is None:
            continue
        loci = copy_loci(h, r.contig, int(r.start), int(r.end), k)
        if not loci:
            continue
        segs.at[i, "copies"] = ",".join(f"{d['contig']}:{d['start']}-{d['end']}:{d['kmers']}" for d in loci[:6])
        if int(r.dominant_m) < 1:
            continue
        ratio, rel = _ratio_of(w, model, r, ok)
        if not np.isfinite(ratio):
            continue
        u = unique_at_copies(w, model, p, ok, _windows_inside(w, index, loci), ratio, rel)
        if u["n"] > 0:
            segs.at[i, "copy_unique"] = f"{u['n']}:{u['copies']}"
            segs.at[i, "copy_llr"] = u["llr"]
        if r.kind != "CN_OR_DEPTH" or u["n"] < min_unique:
            continue
        own = float(r.llr_cn) - float(r.llr_depth)  # the segment's own evidence for a copy change over a factor
        if u["llr"] >= p.margin and u["llr"] + own >= p.margin and abs(u["z_one"]) <= p.max_z:
            segs.at[i, "kind"] = "CN"
            segs.at[i, "resolved"] = "copies"
        elif u["llr"] <= -p.margin and u["llr"] + own <= -p.margin and abs(u["z_ratio"]) <= p.max_z:
            segs.at[i, "kind"] = "DEPTH"
            segs.at[i, "resolved"] = "copies"
    return segs
