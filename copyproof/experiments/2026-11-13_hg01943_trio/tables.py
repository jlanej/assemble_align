#!/usr/bin/env python3
"""The tables of this record, from a finished copyproof run on HG01943.

  tables.py RUN_DIR OUT_DIR

Writes issues_at_locus.tsv, issues_by_class.tsv, joins.tsv (the unlinked joins with what lies next
to them), joins_evidence.tsv (for each join with a copy-number issue beside it: the size of the
gap, the issue, and whether its sequence is repeated elsewhere on the same scaffold) and
missing_in_neither_parent.tsv (one-sided missing sequence whose k-mers are solid in the child and
in neither parent).
"""

import re
import sys
from pathlib import Path

import pandas as pd

run, out = Path(sys.argv[1]), Path(sys.argv[2])
pd.set_option("display.width", 220)
PAT, MAT = "HG01943#1#CM089354.1", "HG01943#2#CM089362.1"

iss = pd.read_csv(run / "issues.tsv", sep="\t")
locus = iss[((iss["contig"] == PAT) & (iss["start"] < 48_700_000) & (iss["end"] > 47_200_000))
            | ((iss["contig"] == MAT) & (iss["start"] < 48_300_000) & (iss["end"] > 46_800_000))]
locus.to_csv(out / "issues_at_locus.tsv", sep="\t", index=False)
iss.groupby(["class", "confidence", "role"]).size().reset_index(name="n").to_csv(out / "issues_by_class.tsv", sep="\t", index=False)
joins = iss[iss["class"] == "UNSUPPORTED_JOIN"]
joins.to_csv(out / "joins.tsv", sep="\t", index=False)
print(f"issues: {len(iss)}; at the locus: {len(locus)}")
print(f"joins without a linking long read: {len(joins)}; next to another issue: {int((joins['confidence'] == 'high').sum())}; "
      f"of those outside the locus: {int(((joins['confidence'] == 'high') & ~joins['id'].isin(locus['id'])).sum())}")
# What is beside each firm join? At the known locus the surplus copy beside a join repeats sequence
# from beside another join of the same scaffold: two pieces that overlap were joined end to end.
gaps = pd.read_csv(run / "scan.gaps.bed", sep="\t", header=None, names=["contig", "start", "end"])
rows = []
for j in joins[joins["confidence"] == "high"].itertuples(index=False):
    g = gaps[(gaps["contig"] == j.contig) & (gaps["start"] < j.end) & (gaps["end"] > j.start)]
    near = iss[(iss["contig"] == j.contig) & (iss["start"] < j.end + 20_000) & (iss["end"] > j.start - 20_000)
               & iss["class"].isin(["EXPANSION", "COLLAPSE"])]
    for n in near.itertuples(index=False):
        ev = dict(kv.split("=", 1) for kv in n.evidence.split(";") if "=" in kv)
        same = []
        for c in ev.get("other_copies", "").split(","):
            m = re.match(r"(.+):(\d+)-(\d+):(\d+)$", c)
            if m and m.group(1) == j.contig:
                same.append(f"{m.group(2)}-{m.group(3)}")
        rows.append({
            "contig": j.contig, "gap_start": int(g["start"].min()) if len(g) else -1,
            "gap_size": int((g["end"] - g["start"]).sum()) if len(g) else -1, "at_known_locus": j.id in set(locus["id"]),
            "issue": n.id, "class": getattr(n, "_4", None) or n[4], "start": int(n.start), "end": int(n.end),
            "confidence": n.confidence, "copies_by_multiplicity": ev.get("copies_by_multiplicity", ""),
            "repeated_on_the_same_sequence": ";".join(same),
        })
je = pd.DataFrame(rows)
je.to_csv(out / "joins_evidence.tsv", sep="\t", index=False)
ex = je[je["class"] == "EXPANSION"]
for known, g in ex.groupby("at_known_locus"):
    rep = g[g["repeated_on_the_same_sequence"] != ""]
    print(f"expansions beside a firm join, {'at the known locus' if known else 'elsewhere'}: {len(g)}; repeated on the same "
          f"sequence: {len(rep)}; gap sizes: {sorted(set(g['gap_size']))}")

s1 = iss[iss["channel"].str.startswith("S1")]
print("ambiguous copy-number segments reported, by support:",
      dict(s1["evidence"].str.extract(r"support=([a-z_]+)")[0].dropna().value_counts()))

cl = pd.read_csv(run / "missing_clusters.tsv", sep="\t")
ms = cl[(cl["kind"] == "MISSING_SEQUENCE") & (cl["frac_in_parent1"] == 0) & (cl["frac_in_parent2"] == 0)].copy()
ms = ms.sort_values("median_count", ascending=False)
cols = ["contig", "start", "n_reads_locus", "reads_fwd", "reads_rev", "tracked_kmers", "median_count", "twin_group", "other_loci"]
ms[cols].to_csv(out / "missing_in_neither_parent.tsv", sep="\t", index=False)
strong = (cl["n_reads_locus"] >= 6) & (cl["tracked_kmers"] >= 5)
allms = cl[cl["kind"] == "MISSING_SEQUENCE"]
print(f"one-sided missing sequence: {len(allms)} placements ({int(strong[allms.index].sum())} strong); in neither parent: "
      f"{len(ms)} placements, {ms['twin_group'].nunique()} events")
print("highest k-mer counts among those in neither parent:")
print(ms[cols[:7]].head(4).to_string(index=False))
