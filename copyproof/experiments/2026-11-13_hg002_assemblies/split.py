#!/usr/bin/env python3
"""High-confidence copy-number issues confirmed by the truth, split by when each assembly was run.

  split.py scoreboard.tsv OUT.tsv

The first seven assemblies were scored while the classifier and the reporting rules were still
being revised (NOTEBOOK entry j), so their scores could shape the rules. The other six were run
after the last change (entry o). All thirteen are scored by the same, final code. One assembly
(`verkko_2023`) holds most of the issues, so each group is also given without it, and the later
group without the two drafts of the truth genome, which have almost nothing left to find.
"""

from __future__ import annotations

import sys

import pandas as pd
from scipy.stats import fisher_exact

FIRST = ["hprc_year2", "verkko_2023", "hprc_year1", "hifiasm_2021", "hifiasm_ontonly_2025", "hg002v0.7", "hprc_year2_v2_thic"]
BIG = "verkko_2023"
DRAFTS = ["hg002v0.9", "hg002v1.0.1"]

t = pd.read_csv(sys.argv[1], sep="\t").set_index("assembly")
later = [a for a in t.index if a not in FIRST]
groups = [
    ("scored while the rules were still being revised", FIRST),
    (f"... without {BIG}", [a for a in FIRST if a != BIG]),
    ("run after the last change", later),
    ("... without the two drafts", [a for a in later if a not in DRAFTS]),
    (f"all without {BIG}", [a for a in t.index if a != BIG]),
    ("all", list(t.index)),
]
rows = []
for label, names in groups:
    high, true = int(t.loc[names, "high"].sum()), int(t.loc[names, "high_true"].sum())
    rows.append({"group": label, "assemblies": len(names), "high": high, "confirmed": true,
                 "not_confirmed": high - true, "share": round(true / high, 4)})
out = pd.DataFrame(rows)
out.to_csv(sys.argv[2], sep="\t", index=False)
print(out.to_string(index=False))
a, b = rows[1], rows[2]
p = fisher_exact([[a["confirmed"], a["not_confirmed"]], [b["confirmed"], b["not_confirmed"]]])[1]
print(f"Fisher's exact test, '{a['group']}' against '{b['group']}': p = {p:.2f}")
