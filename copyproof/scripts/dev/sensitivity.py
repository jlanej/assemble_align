#!/usr/bin/env python3
"""How much do the results depend on the classifier's thresholds, and on its revision?

  sensitivity.py --truth RUN --planted RUN --expected plant.expected.bed \
                 --asm RUN --asm-truth NAME.eval.truth_cn_segments_raw.tsv --out OUT.tsv

Each run directory must hold windows, a model and a cached segmentation (scripts/dev/cache_segments.py).
The HMM is not run again: only the classification of its segments and the rules that turn
segments into issues. For every setting the same three things are measured:

  truth genome    whole-copy segments, and copy-number issues of high and of low confidence
                  (on a truth genome fewer is better)
  planted genome  planted duplications, extra copies and deletions found as a copy-number issue,
                  at any confidence and at high confidence
  one assembly    high-confidence issues confirmed by the truth genome's k-mers, and assessable
                  discordances found

Two steps of the pipeline are left out, so the numbers are a little below the pipeline's: the
lookup of a segment's other copies, and support from the other channels.
"""

from __future__ import annotations

import argparse
import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from classifier_before_revision import classify_segment_before  # noqa: E402

from copyproof import calls as C  # noqa: E402
from copyproof.cn import CnParams, assessable, classify_segment  # noqa: E402
from copyproof.model import fit_model  # noqa: E402
from copyproof.truth import match_segments  # noqa: E402
from copyproof.windows import Windows  # noqa: E402

# (label, classifier, parameter overrides). The first two rows are the revision itself.
SETTINGS = [("before the revision", "before", {}), ("as shipped", "final", {})]
for key, values in (
    ("psi_local", (0.03, 0.06)), ("max_z", (3.0, 4.0)), ("margin", (2.0, 4.0)), ("dominant_tol", (0.25, 0.45)),
    ("fit_share", (0.8, 1.0)), ("max_wander", (1.5, 3.0)), ("min_share", (0.4, 0.6)), ("min_seg_class", (10, 40)),
    ("min_llr", (30.0, 40.0)), ("deep_windows", (5, 12)), ("deep_ratio", (1.33, 2.0)), ("high_windows", (3, 8)),
):
    for v in values:
        SETTINGS.append((f"{key} = {v:g}", "final", {key: v}))


class Run:
    def __init__(self, path):
        self.path = Path(path)
        self.w = Windows.load(str(sorted(self.path.glob("windows.w*.npz"))[0]))
        self.model = fit_model(self.w)
        self.ok = assessable(self.w, self.model, CnParams())
        self.cache = pickle.load(open(self.path / "segcache.pkl", "rb"))

    def segments(self, which: str, p: CnParams) -> pd.DataFrame:
        classify = classify_segment_before if which == "before" else classify_segment
        w, rows = self.w, []
        for ci, a, b, d, llr, parts in self.cache:
            if llr < p.min_llr:
                continue
            sl = w.slice_of(ci)
            name, L = w.contigs.at[ci, "name"], int(w.contigs.at[ci, "length"])

            def add(a, b, d, llr):
                w0, w1 = sl.start + a, sl.start + b
                row = {"contig": name, "start": a * w.W, "end": min(b * w.W, L), "direction": "gain" if d > 0 else "loss",
                       "llr": round(llr, 2), "n_win": b - a, "n_assessed": int(self.ok[w0:w1].sum())}
                row.update(classify(w, self.model, p, w0, w1, self.ok))
                return row

            row = add(a, b, d, llr)
            if row["kind"] == "IRREGULAR":
                ps = [add(x, y, d, l) for x, y, l in parts if l >= p.min_llr]
                if len(ps) > 1 and [r for r in ps if r["kind"] != "IRREGULAR"]:
                    rows.extend(ps)
                    continue
            rows.append(row)
        segs = pd.DataFrame(rows)
        if "own_classes" not in segs:
            segs["own_classes"] = ""
        return segs

    def issues(self, segs: pd.DataFrame, p: CnParams) -> pd.DataFrame:
        iss = C.issues_from_cn(segs, False, None, p.deep_ratio, p.deep_windows, p.high_windows)
        iss["direction"] = np.where(iss["class"] == "COLLAPSE", "gain", "loss")
        return iss


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--planted", required=True)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--asm", required=True)
    ap.add_argument("--asm-truth", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    truth, planted, asm = Run(a.truth), Run(a.planted), Run(a.asm)
    exp = pd.read_csv(a.expected, sep="\t", names=["contig", "start", "end", "signal", "tag"])
    exp = exp[exp["signal"].isin(["EXPANSION", "EXPANSION_PARTNER", "COLLAPSE"])]
    tseg = pd.read_csv(a.asm_truth, sep="\t")
    mask = pd.read_csv(asm.path / "cannot_assess.bed", sep="\t", names=["contig", "start", "end"])
    mby = {c: g[["start", "end"]].values for c, g in mask.groupby("contig")}

    def masked(c, s, e):
        m = mby.get(c)
        if m is None:
            return 0.0
        ov = np.minimum(m[:, 1], e) - np.maximum(m[:, 0], s)
        return float(ov[ov > 0].sum()) / max(e - s, 1)

    tseg["masked"] = [masked(c, s, e) for c, s, e in zip(tseg["contig"], tseg["start"], tseg["end"])]

    rows = []
    for label, which, over in SETTINGS:
        p = dataclasses.replace(CnParams(), **over)
        row = {"setting": label}
        # truth genome
        segs = truth.segments(which, p)
        iss = truth.issues(segs, p)
        row["truth_whole_copy_segments"] = int((segs["kind"] == "CN").sum())
        row["truth_ambiguous_segments"] = int((segs["kind"] == "CN_OR_DEPTH").sum())
        row["truth_issues_high"] = int((iss["confidence"] == "high").sum())
        row["truth_issues_low"] = int((iss["confidence"] == "low").sum())
        # planted genome
        segs = planted.segments(which, p)
        iss = planted.issues(segs, p)
        by = {c: g for c, g in iss.groupby("contig")}
        found = high = 0
        for r in exp.itertuples(index=False):
            g = by.get(r.contig)
            if g is None:
                continue
            want = ["COLLAPSE"] if r.signal == "COLLAPSE" else ["EXPANSION", "UNSUPPORTED"]
            o = g[g["class"].isin(want) & (g["start"] < r.end) & (g["end"] > r.start)]
            cov = float((np.minimum(o["end"], r.end) - np.maximum(o["start"], r.start)).clip(lower=0).sum()) / max(r.end - r.start, 1)
            found += cov >= 0.5
            high += cov >= 0.5 and bool((o["confidence"] == "high").any())
        row["planted_signals"] = len(exp)
        row["planted_found"] = int(found)
        row["planted_found_high"] = int(high)
        # one assembly against the truth
        segs = asm.segments(which, p)
        iss = asm.issues(segs, p)
        iss, t = match_segments(iss, tseg)
        ok = t[t["masked"] < 0.5]
        hi = iss[iss["confidence"] == "high"]
        row["asm_high"] = len(hi)
        row["asm_high_confirmed"] = int(hi["true"].sum())
        row["asm_low"] = int((iss["confidence"] == "low").sum())
        row["asm_low_confirmed"] = int(iss.loc[iss["confidence"] == "low", "true"].sum())
        for d, name in (("gain", "collapses"), ("loss", "surplus")):
            row[f"asm_{name}"] = int((ok["direction"] == d).sum())
            row[f"asm_{name}_found"] = int(ok.loc[ok["direction"] == d, "found"].sum())
        rows.append(row)
        print(row, flush=True)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
