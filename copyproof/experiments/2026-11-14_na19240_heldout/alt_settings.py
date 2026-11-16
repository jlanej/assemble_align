#!/usr/bin/env python3
"""The two candidate values, scored on a finished run (NOTEBOOK entry p, "the two candidate values").

`min_seg_class` 10 and `margin` 2 looked better than the shipped 20 and 3 on the data the method was
tuned on (../2026-11-13_threshold_sensitivity/). Here the segments of a run are classified again
under the four combinations, as that record did: the HMM is not run again, and the lookup of other
copies and the support from other channels are left out, so "as shipped" sits a little below the
pipeline's own numbers. Every issue list is scored against the discordance tables that
scripts/eval_assemblies.py wrote for the run (PREFIX.segments.tsv), one row per setting and label set.

  ../../scripts/dev/cache_segments.py RUN_DIR          # once: the HMM's segments
  alt_settings.py RUN_DIR EVAL_PREFIX OUT.tsv
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve()
sys.path.insert(0, str(here.parents[2] / "src"))
sys.path.insert(0, str(here.parents[2] / "scripts" / "dev"))
from sensitivity import Run  # noqa: E402

from copyproof.cn import CnParams  # noqa: E402
from copyproof.truth import match_segments  # noqa: E402

SETTINGS = [("as shipped (20, 3)", {}), ("min_seg_class = 10", {"min_seg_class": 10}), ("margin = 2", {"margin": 2.0}),
            ("min_seg_class = 10 and margin = 2", {"min_seg_class": 10, "margin": 2.0})]


def main() -> None:
    run_dir, prefix, out = sys.argv[1:4]
    run = Run(run_dir)
    seg = pd.read_csv(f"{prefix}.segments.tsv", sep="\t")
    rows = []
    for label, over in SETTINGS:
        p = dataclasses.replace(CnParams(), **over)
        segs = run.segments("final", p)
        iss = run.issues(segs, p)
        for name, t in seg.groupby("labels", sort=False):
            t = t.drop(columns=["call_overlap", "found", "found_high"])
            scored, tt = match_segments(iss, t)
            hi, lo = scored[scored["confidence"] == "high"], scored[scored["confidence"] == "low"]
            ok = tt[(tt["masked"] < 0.5) & (tt["len"] >= 5000)]
            row = {"setting": label, "labels": name, "whole_copy_segments": int((segs["kind"] == "CN").sum()),
                   "ambiguous_segments": int((segs["kind"] == "CN_OR_DEPTH").sum()), "high": len(hi),
                   "high_confirmed": int(hi["true"].sum()), "low": len(lo), "low_confirmed": int(lo["true"].sum())}
            for d, what in (("gain", "collapses"), ("loss", "surplus")):
                row[f"{what}_5kb"] = int((ok["direction"] == d).sum())
                row[f"{what}_5kb_found"] = int(ok.loc[ok["direction"] == d, "found"].sum())
            rows.append(row)
            print(row, flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
