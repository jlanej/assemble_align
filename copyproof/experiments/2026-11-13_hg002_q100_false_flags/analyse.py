#!/usr/bin/env python3
"""Rung 3: what does copyproof flag on the truth genome itself?

  analyse.py RUN_DIR ANNOT_DIR OUT_PREFIX

RUN_DIR is a copyproof run on hg002v1.1 with HG002's own reads. ANNOT_DIR holds the Q100
annotation (issues v4, SDs, centromeres, HSats). Writes OUT_PREFIX.*.tsv and prints the tables.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402

run, annot, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 400)
pd.set_option("display.max_colwidth", 90)

summary = json.loads((run / "summary.json").read_text())
cen = IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.centromeres.bed"))
hsat = IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.hsats.bed"))
sd = IntervalSet.from_bed(str(annot / "v4_issue_breakdown/hg002v1.1.SDs.013025.merged.bed"))
known = IntervalSet.from_bed(str(annot / "hg002v1.1_issues.v4.bed"))
censat = cen.union(hsat)


def stratum(c, s, e):
    L = max(e - s, 1)
    if censat.overlap(c, s, e) > 0.5 * L:
        return "censat"
    if sd.overlap(c, s, e) > 0.5 * L:
        return "SD"
    return "unique"


print("== run")
for k in ("depth_per_copy", "assessable_fraction", "cannot_assess_bp", "kmer_support", "read_kmer_completeness", "cn_segments",
          "unsupported_clusters", "missing", "issues", "timings"):
    print(f"{k}: {summary.get(k)}")
print("depth model:", {k: summary["depth_model"][k] for k in ("omega", "psi", "ratio_sd_per_window")})

# the mask by stratum
mask = pd.read_csv(run / "cannot_assess.bed", sep="\t", names=["contig", "start", "end"])
contigs = pd.read_csv(run / "scan.contigs.tsv", sep="\t")
total = {"censat": censat.total(), "SD": sd.total()}
m = IntervalSet.from_rows(zip(mask.contig, mask.start, mask.end))
masked = {"all": m.total()}
for name, st in (("censat", censat), ("SD", sd)):
    masked[name] = sum(m.overlap(c, s, e) for c, s, e in st.to_rows())
genome = int(contigs["length"].sum())
print("\n== cannot assess")
print(f"genome {genome / 1e6:.1f} Mb; masked {masked['all'] / 1e6:.1f} Mb ({100 * masked['all'] / genome:.2f}%)")
for name in ("censat", "SD"):
    print(f"  of {name} ({total[name] / 1e6:.1f} Mb): {100 * masked[name] / total[name]:.1f}% masked")

segs = pd.read_csv(run / "cn_segments.tsv", sep="\t")
segs["len"] = segs["end"] - segs["start"]
segs["stratum"] = [stratum(c, s, e) for c, s, e in zip(segs.contig, segs.start, segs.end)]
segs["known_bp"] = [known.overlap(c, s, e) for c, s, e in zip(segs.contig, segs.start, segs.end)]
print("\n== copy-number segments by kind")
print(segs.groupby(["kind", "direction", "stratum"]).agg(n=("len", "size"), Mb=("len", lambda x: round(x.sum() / 1e6, 2)),
                                                           median_len=("len", "median"), known=("known_bp", lambda x: int((x > 0).sum()))).to_string())
segs.to_csv(f"{out}.cn_segments_annotated.tsv", sep="\t", index=False)
anom = segs[segs["kind"].isin(["DEPTH", "IRREGULAR"])]
print(f"depth anomalies: {len(anom)}; on a known issue: {int((anom['known_bp'] > 0).sum())}")
cand = segs[segs["kind"] == "CN"].sort_values("llr", ascending=False)
cand[["contig", "start", "end", "direction", "delta", "own", "llr", "stratum", "known_bp", "classes", "own_classes"]].to_csv(
    f"{out}.cn_candidates.tsv", sep="\t", index=False)

iss = pd.read_csv(run / "issues.tsv", sep="\t")
iss["len"] = iss["end"] - iss["start"]
iss["stratum"] = [stratum(c, s, e) for c, s, e in zip(iss.contig, iss.start, iss.end)]
iss["known_bp"] = [known.overlap(c, max(s - 100, 0), e + 100) for c, s, e in zip(iss.contig, iss.start, iss.end)]
iss["group"] = np.where(iss.contig.str.contains("MATERNAL"), "mat", np.where(iss.contig.str.contains("PATERNAL"), "pat", "other"))
print("\n== issues on the truth genome")
t = iss.groupby(["class", "confidence", "role"]).agg(n=("len", "size"), in_known_issue=("known_bp", lambda x: int((x > 0).sum())),
                                                     unique=("stratum", lambda x: int((x == "unique").sum())),
                                                     SD=("stratum", lambda x: int((x == "SD").sum())),
                                                     censat=("stratum", lambda x: int((x == "censat").sum())))
print(t.to_string())
t.to_csv(f"{out}.issues_by_class.tsv", sep="\t")
iss.to_csv(f"{out}.issues_annotated.tsv", sep="\t", index=False)
print("\nper haplotype (all classes / high confidence, role locus or unresolved):")
core = iss[iss["role"].isin(["locus", "unresolved"])]
print(core.groupby("group").agg(all=("len", "size"), high=("confidence", lambda x: int((x == "high").sum()))).to_string())

cols = ["id", "contig", "start", "len", "class", "role", "confidence", "score", "stratum", "known_bp", "evidence"]
for cls in ("COLLAPSE", "EXPANSION", "UNSUPPORTED", "MISSING_JUNCTION", "MISSING_HAPLOTYPE", "FALSE_JOIN", "MISSING_SEQUENCE"):
    sub = iss[iss["class"] == cls].sort_values(["confidence", "score"], ascending=[True, False])
    if len(sub):
        print(f"\n== {cls}: {len(sub)}")
        print(sub.head(12)[cols].to_string(index=False))
