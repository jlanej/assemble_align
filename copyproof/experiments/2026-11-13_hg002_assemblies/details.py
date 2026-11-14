#!/usr/bin/env python3
"""Rung 4, two questions the scoreboard does not answer, for one assembly.

  details.py WORK NAME OUT_PREFIX

WORK holds NAME/ (the copyproof run) and NAME.eval.truth_cn_segments.tsv (from eval_truth.py).

1. What does copyproof see at the collapses it misses? For every assessable discordance where the
   truth holds more copies than the assembly: the segments that overlap it, by kind.
2. Are short whole-copy calls real? Whole-copy segments on fewer than 5 assessed windows, scored
   against the truth whatever confidence the issue list gives them.
3. If the run had the reads: are the clusters of read sequence that the assembly lacks real? A
   cluster is true when the truth genome holds at least half of its missing k-mers; strong means
   at least 6 reads at the locus and 5 tracked k-mers, as in the issue list.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.truth import match_segments  # noqa: E402

work, name, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
pd.set_option("display.width", 220)
truth = pd.read_csv(work / f"{name}.eval.truth_cn_segments.tsv", sep="\t")
segs = pd.read_csv(work / name / "cn_segments.tsv", sep="\t")
by = {c: g for c, g in segs.groupby("contig")}

rows = []
for r in truth[(truth["direction"] == "gain") & (truth["masked"] < 0.5)].itertuples(index=False):
    g = by.get(r.contig)
    ov = g[(g["start"] < r.end) & (g["end"] > r.start)] if g is not None else segs.iloc[:0]
    cover = {}
    for kind, x in ov.groupby("kind"):
        cover[kind] = float((np.minimum(x["end"], r.end) - np.maximum(x["start"], r.start)).sum()) / (r.end - r.start)
    best = max(cover, key=cover.get) if cover else "none"
    top = ov.assign(o=np.minimum(ov["end"], r.end) - np.maximum(ov["start"], r.start)).sort_values("o", ascending=False).head(1)
    rows.append({"contig": r.contig, "start": r.start, "end": r.end, "found": bool(r.found), "mean_delta": r.mean_delta,
                 "main_segment_kind": best, "covered": round(sum(cover.values()), 2),
                 "classes": top.iloc[0]["classes"] if len(top) else ""})
col = pd.DataFrame(rows)
col.to_csv(f"{out}.collapses.tsv", sep="\t", index=False)
print(f"{name}: assessable discordances where the truth holds more copies: {len(col)}; found: {int(col['found'].sum())}")
print("the others, by the kind of segment that covers most of them:")
print(col[~col["found"]].groupby("main_segment_kind").size().to_string())

cn = segs[segs["kind"] == "CN"].copy()
cn, _ = match_segments(cn, truth)
cn["windows"] = np.where(cn["n_assessed"] < 5, "under 5", "5 or more")
t = cn.groupby(["direction", "windows"]).agg(segments=("true", "size"), confirmed=("true", "sum")).reset_index()
t.to_csv(f"{out}.whole_copy_by_length.tsv", sep="\t", index=False)
print("\nwhole-copy segments against the truth, by assessed windows:")
print(t.to_string(index=False))
short = cn[(cn["windows"] == "under 5") & (cn["direction"] == "loss")]
print(f"short losses: {len(short)}, confirmed {int(short['true'].sum())}; on the haplotype's own k-mers: "
      f"{int((short['own'] == 'own').sum())}, confirmed {int(short.loc[short['own'] == 'own', 'true'].sum())}")

mc = work / name / "missing_clusters.tsv"
if mc.is_file():
    c = pd.read_csv(mc, sep="\t")
    if "frac_in_truth" in c.columns and len(c):
        c["strength"] = np.where((c["n_reads_locus"] >= 6) & (c["tracked_kmers"] >= 5), "strong", "weak")
        c["true"] = c["frac_in_truth"] >= 0.5
        m = c.groupby(["kind", "strength"]).agg(clusters=("true", "size"), true=("true", "sum")).reset_index()
        m["share_true"] = (m["true"] / m["clusters"]).round(3)
        m.to_csv(f"{out}.missing_clusters.tsv", sep="\t", index=False)
        print("\nclusters of read sequence the assembly lacks, against the truth genome's k-mers:")
        print(m.to_string(index=False))
