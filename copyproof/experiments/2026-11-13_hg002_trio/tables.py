#!/usr/bin/env python3
"""The hapmer channel on a truth genome and on planted haplotype swaps.

  tables.py TRUTH_TRIO_RUN PLANTED_TRIO_RUN plant.expected.bed ANNOT_DIR OUT_DIR

Both runs are copyproof runs with the parents' read sets (scripts/run_hg002_trio.sh). The audit is
repeated here with and without its two filters, so the effect of each is measured on the same
windows: on the truth genome every switch is a false flag, on the planted genome the swapped
stretches are the positives.
"""

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.hapmers import HapParams, audit, readable_windows  # noqa: E402
from copyproof.intervals import IntervalSet  # noqa: E402
from copyproof.track import Scan  # noqa: E402
from copyproof.windows import Windows  # noqa: E402

truth_run, plant_run, expected, annot, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), Path(sys.argv[5])
pd.set_option("display.width", 220)
OFF = dict(parent_depth=(0.0, 1e9), min_density=0.0)
SETTINGS = (("no filter", OFF), ("parents readable", dict(min_density=0.0)), ("site density", dict(parent_depth=(0.0, 1e9))),
            ("both (as shipped)", {}))


def load(run):
    sc = Scan(str(run / "scan"))
    w = Windows.load(str(sorted(run.glob("windows.w*.npz"))[0]))
    return sc, w


def switches(sc, w, over):
    p = dataclasses.replace(HapParams(), **over)
    contigs, pieces, sw = audit(w, sc.gaps, (sc.reads[1], sc.reads[2]), p)
    if len(sw):
        sw["high"] = sw["sites_segment_parent"] >= 10
    return contigs, pieces, sw


# ---- truth genome: every switch is a false flag ----------------------------------------------------
sc, w = load(truth_run)
p = HapParams()
ok = readable_windows(w, p)
print(f"truth genome: windows where a parent's depth is outside {p.parent_depth} of its usual depth: "
      f"{int((~ok).sum())} of {len(ok)} ({100 * (~ok).mean():.2f}%)")
rows = []
for label, over in SETTINGS:
    contigs, pieces, sw = switches(sc, w, over)
    rows.append({"setting": label, "switch_segments": len(sw), "high_confidence": int(sw["high"].sum()) if len(sw) else 0,
                 "Mb": round(float((sw["end"] - sw["start"]).sum()) / 1e6, 2) if len(sw) else 0.0,
                 "pieces_other_parent": int((pieces["call"] == "other_parent").sum()) if len(pieces) else 0})
truth_table = pd.DataFrame(rows)
print(truth_table.to_string(index=False))
contigs, pieces, sw = switches(sc, w, {})
assessed = contigs[contigs["parent"] != "not_assessable"]
print(f"sequences assessed: {len(assessed)}; hapmer sites: {int(assessed[[c for c in contigs.columns if c.startswith('sites_') and c != 'sites_neither']].sum().sum())}; "
      f"purity {assessed['purity'].min():.4f} to {assessed['purity'].max():.4f}; pieces: {len(pieces)}, "
      f"{dict(pieces['call'].value_counts()) if len(pieces) else {}}")
cen = IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.centromeres.bed")).union(
    IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.hsats.bed")))
sd = IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.SDs.013025.merged.bed"))
known = IntervalSet.from_bed(str(annot / "hg002v1.1_issues.v4.bed"))
if len(sw):
    L = (sw["end"] - sw["start"]).clip(lower=1)
    sw["stratum"] = np.where([cen.overlap(c, s, e) for c, s, e in zip(sw.contig, sw.start, sw.end)] > 0.5 * L, "censat",
                             np.where([sd.overlap(c, s, e) for c, s, e in zip(sw.contig, sw.start, sw.end)] > 0.5 * L, "SD", "unique"))
    sw["known_issue"] = [known.overlap(c, s, e) > 0 for c, s, e in zip(sw.contig, sw.start, sw.end)]
    sw.sort_values("sites_segment_parent", ascending=False).to_csv(out / "truth_switches.tsv", sep="\t", index=False)
    print("\nswitch segments kept on the truth genome:")
    print(sw.groupby(["stratum", "high"]).agg(n=("llr", "size"), known_issue=("known_issue", "sum"),
                                              median_kb=("start", lambda x: float(np.median((sw.loc[x.index, "end"] - x) / 1000)))).to_string())
    print(sw.sort_values("sites_segment_parent", ascending=False).head(10)[
        ["contig", "start", "end", "segment_parent", "sites_segment_parent", "sites_scaffold_parent", "stratum", "known_issue"]].to_string(index=False))
truth_table.to_csv(out / "truth_filters.tsv", sep="\t", index=False)

# what the parents' reads say about read sequence the truth genome lacks
cl = pd.read_csv(truth_run / "missing_clusters.tsv", sep="\t")
ms = cl[cl["kind"] == "MISSING_SEQUENCE"]
strong = (ms["n_reads_locus"] >= 6) & (ms["tracked_kmers"] >= 5)
inherited = ms[["frac_in_parent1", "frac_in_parent2"]].max(axis=1) >= 0.5
neither = (ms["frac_in_parent1"] == 0) & (ms["frac_in_parent2"] == 0)
print(f"\none-sided missing sequence on the truth genome: {len(ms)} placements; strong {int(strong.sum())}, of which a parent's "
      f"reads hold most k-mers: {int((strong & inherited).sum())}; in neither parent: {int(neither.sum())} placements "
      f"({ms.loc[neither, 'twin_group'].nunique()} events), of which weak: {int((neither & ~strong).sum())}")

# ---- planted genome: the swapped stretches -----------------------------------------------------------
sc2, w2 = load(plant_run)
exp = pd.read_csv(expected, sep="\t", names=["contig", "start", "end", "signal", "tag"])
exp = exp[exp["signal"] == "PHASE_SWITCH"].copy()
exp["size"] = exp["tag"].str.split("|").str[2].astype(int)
rows = []
for label, over in SETTINGS:
    _, _, s2 = switches(sc2, w2, over)
    by = {c: g for c, g in s2.groupby("contig")} if len(s2) else {}
    for size, g in exp.groupby("size"):
        found = 0
        for r in g.itertuples(index=False):
            x = by.get(r.contig)
            if x is None:
                continue
            o = x[(x["start"] < r.end) & (x["end"] > r.start)]
            cov = float((np.minimum(o["end"], r.end) - np.maximum(o["start"], r.start)).clip(lower=0).sum()) / (r.end - r.start)
            found += cov >= 0.5
        rows.append({"setting": label, "swap_size": size, "planted": len(g), "found": int(found), "switch_segments_in_all": len(s2)})
plant_table = pd.DataFrame(rows)
plant_table.to_csv(out / "planted_swaps.tsv", sep="\t", index=False)
print("\nplanted haplotype swaps found as a phase switch (at least half covered):")
print(plant_table.pivot(index="setting", columns="swap_size", values="found").reindex([s for s, _ in SETTINGS]).to_string())
print(plant_table.groupby("setting")["switch_segments_in_all"].first().reindex([s for s, _ in SETTINGS]).to_string())
