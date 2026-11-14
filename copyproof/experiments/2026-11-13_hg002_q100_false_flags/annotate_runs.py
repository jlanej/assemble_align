#!/usr/bin/env python3
"""Rung 3 with long reads: are the runs that neither read type holds known problems of the truth genome?

  annotate_runs.py RUN_DIR ANNOT_DIR OUT.tsv

RUN_DIR is a copyproof run on hg002v1.1 with --long-read-bam (it holds unsupported_points.tsv and
unsupported_clusters.tsv). ANNOT_DIR holds the Q100 annotation: hg002v1.1_issues.v4.bed and
changes/hg002v1.1_to_hg002v1.2.vcf.gz. Every audited run of unsupported assembly k-mers gets two
flags: on a known-issue interval (within 100 bp), and within 31 bp of a base that v1.2 changed.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

run, annot, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
pd.set_option("display.width", 220)

points = pd.read_csv(run / "unsupported_points.tsv", sep="\t", dtype={"name": str})
runs = points[points["name"].str.startswith("run")].copy()
losses = points[points["name"].str.startswith("cn")]

rows = []
with gzip.open(annot / "changes/hg002v1.1_to_hg002v1.2.vcf.gz", "rt") as f:
    for line in f:
        if not line.startswith("#"):
            t = line.split("\t", 5)
            rows.append((t[0], int(t[1]) - 1, int(t[1]) - 1 + len(t[3])))
changed = pd.DataFrame(rows, columns=["contig", "start", "end"])
known = pd.read_csv(annot / "hg002v1.1_issues.v4.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["contig", "start", "end"])


def index(t):
    return {c: g[["start", "end"]].values for c, g in t.groupby("contig")}


def near(tab, c, s, e, pad):
    r = tab.get(c)
    return r is not None and bool(np.any((r[:, 0] < e + pad) & (r[:, 1] > s - pad)))


cb, kb = index(changed), index(known)
runs["fixed_in_v1.2"] = [near(cb, c, s, e, 31) for c, s, e in zip(runs["contig"], runs["start"], runs["end"])]
runs["known_issue"] = [near(kb, c, s, e, 100) for c, s, e in zip(runs["contig"], runs["start"], runs["end"])]
runs.to_csv(out, sep="\t", index=False)

print(f"changes from v1.1 to v1.2: {len(changed)} ({int((changed['end'] - changed['start']).sum())} bases); "
      f"known-issue intervals: {len(known)}")
print("\nruns of unsupported assembly k-mers with covered flanks, by what long reads hold:")
t = runs.groupby("long_read_kmers").agg(runs=("name", "size"), known_issue=("known_issue", "sum"), changed_in_v12=("fixed_in_v1.2", "sum"))
print(t.to_string())
nc = runs[runs["long_read_kmers"] == "not_carried"]
print("\nnot carried, by the verdict of alignment alone:")
print(nc.groupby("verdict").agg(runs=("name", "size"), known_issue=("known_issue", "sum"), changed_in_v12=("fixed_in_v1.2", "sum"),
                                median_spanning_reads=("spanning", "median")).to_string())
new = nc[~nc["fixed_in_v1.2"] & ~nc["known_issue"]]
print(f"\nnot carried, on no known issue and unchanged in v1.2: {len(new)}")
print(new[["contig", "start", "end", "spanning", "clipped", "kmers", "kmers_carried", "long_reads_over"]].to_string(index=False))
print(f"\ncopy-number losses on a haplotype's own k-mers: {dict(losses['long_read_kmers'].value_counts())}")

# how many of the v1.2 corrections does this channel see at all?
clusters = pd.read_csv(run / "unsupported_clusters.tsv", sep="\t")
ub, rb, nb = index(clusters), index(runs), index(nc)
for label, tab in (("any run of unsupported k-mers", ub), ("an audited run (covered flanks)", rb), ("a run long reads do not carry", nb)):
    n = sum(near(tab, c, s, e, 31) for c, s, e in zip(changed["contig"], changed["start"], changed["end"]))
    print(f"v1.2 changes within 31 bp of {label}: {n} of {len(changed)}")
