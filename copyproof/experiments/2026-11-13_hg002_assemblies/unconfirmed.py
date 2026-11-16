#!/usr/bin/env python3
"""Every high-confidence copy-number issue that the truth genome does not confirm, over all assemblies.

  unconfirmed.py WORK assemblies.tsv OUT.tsv

WORK holds NAME.eval.cn_issues_scored.tsv from scripts/eval_truth.py for the names in
assemblies.tsv. `truth_overlap` is the share of the issue covered by a truth discordance of the
same direction (confirmed from 0.5).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

work, names, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
names = [line.split("\t")[0] for line in Path(names).read_text().splitlines() if line and not line.startswith("#")]
rows = []
for name in sorted(names):
    f = work / f"{name}.eval.cn_issues_scored.tsv"
    if not f.is_file():
        continue
    t = pd.read_csv(f, sep="\t")
    u = t[(t["confidence"] == "high") & ~t["true"].astype(bool)].copy()
    u.insert(0, "assembly", name)
    u["kb"] = ((u["end"] - u["start"]) / 1000).round().astype(int)
    rows.append(u[["assembly", "contig", "start", "end", "kb", "class", "truth_overlap", "summary"]])
u = pd.concat(rows, ignore_index=True)
u["truth_overlap"] = u["truth_overlap"].round(2)
u.to_csv(out, sep="\t", index=False)
n = u.groupby("assembly").size()
print(f"{len(u)} unconfirmed high-confidence copy-number issues over {len(rows)} assemblies")
print(n.to_string())
print(f"partly covered by a truth discordance (0.2 to 0.5): {int(((u.truth_overlap >= 0.2) & (u.truth_overlap < 0.5)).sum())}; "
      f"not at all (0): {int((u.truth_overlap == 0).sum())}")
