#!/usr/bin/env python3
"""The false-flag litmus on a truth genome (GUIDE section 8, rung 3; decision D5).

Denominator, fixed before the run: every boundary of an annotated segmental duplication whose
flank outside the duplication is unique for at least 5 kb (no other duplication, no satellite or
centromere annotation, not a sequence end). A boundary is flagged if an issue lies within 5 kb.
The report's target is at most 2% of boundaries flagged.

  sd_boundary_litmus.py --sd SD.bed --avoid censat.bed [more.bed] --fai asm.fai --issues issues.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sd", required=True)
ap.add_argument("--avoid", nargs="*", default=[])
ap.add_argument("--fai", required=True)
ap.add_argument("--issues", required=True)
ap.add_argument("--flank", type=int, default=5000)
ap.add_argument("--near", type=int, default=5000)
ap.add_argument("--out")
a = ap.parse_args()

sd = IntervalSet.from_bed(a.sd)
avoid = IntervalSet()
for p in a.avoid:
    avoid = avoid.union(IntervalSet.from_bed(p))
lengths = {t[0]: int(t[1]) for t in (line.split("\t") for line in open(a.fai))}

bounds = []
for contig, iv in sd.by.items():
    L = lengths.get(contig)
    if L is None:
        continue
    for i, (s, e) in enumerate(iv):
        # left boundary: unique flank is [s - flank, s)
        prev_end = iv[i - 1][1] if i > 0 else 0
        if s - prev_end >= a.flank and s >= a.flank and avoid.overlap(contig, s - a.flank, s + a.flank) == 0:
            bounds.append((contig, int(s), "left"))
        nxt = iv[i + 1][0] if i + 1 < len(iv) else L
        if nxt - e >= a.flank and L - e >= a.flank and avoid.overlap(contig, e - a.flank, e + a.flank) == 0:
            bounds.append((contig, int(e), "right"))
b = pd.DataFrame(bounds, columns=["contig", "pos", "side"])

iss = pd.read_csv(a.issues, sep="\t")
classes = ["COLLAPSE", "EXPANSION", "UNSUPPORTED", "MISSING_JUNCTION", "FALSE_JOIN", "MISSING_HAPLOTYPE", "PHASE_SWITCH",
           "MISPLACED_PIECE", "UNSUPPORTED_JOIN"]
iss = iss[iss["class"].isin(classes)]
by = {c: g for c, g in iss.groupby("contig")}


def flagged(contig, pos, conf=None):
    g = by.get(contig)
    if g is None:
        return ""
    if conf:
        g = g[g["confidence"] == conf]
    hit = g[(g["start"] <= pos + a.near) & (g["end"] >= pos - a.near)]
    return ",".join(sorted(set(hit["class"])))


b["issues"] = [flagged(c, p) for c, p in zip(b["contig"], b["pos"])]
b["issues_high"] = [flagged(c, p, "high") for c, p in zip(b["contig"], b["pos"])]
n = len(b)
k_any, k_high = int((b["issues"] != "").sum()), int((b["issues_high"] != "").sum())
print(f"anchored duplication boundaries: {n}")
print(f"flagged by any issue within {a.near} bp: {k_any} ({100 * k_any / max(n, 1):.2f}%)")
print(f"flagged by a high-confidence issue: {k_high} ({100 * k_high / max(n, 1):.2f}%)")
print(b.loc[b["issues"] != "", "issues"].value_counts().head(12).to_string())
if a.out:
    b.to_csv(a.out, sep="\t", index=False)
