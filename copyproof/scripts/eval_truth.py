#!/usr/bin/env python3
"""Score a copyproof run against a truth genome's k-mers (GUIDE section 8, rungs 3 and 4).

The run must have been scanned with the truth genome's exact k-mer counts as its last read set
(see copyproof/truth.py). Every tracked k-mer of the assembly then carries the number of copies the
truth genome holds, and labels need no alignment:

  copy number   windows where most k-mers have more (or fewer) copies in the truth than in the
                assembly are true discordances; copyproof's CN segments are scored against them
  missing       a solid read k-mer absent from the assembly is a true missing k-mer if the truth
                genome holds it
  unsupported   an assembly k-mer without read support is a true error if the truth lacks it

  eval_truth.py --run RUN_DIR --out PREFIX [--truth-label truth]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from copyproof.model import Model  # noqa: E402
from copyproof.track import Scan  # noqa: E402
from copyproof.truth import match_segments, truth_segments, truth_windows  # noqa: E402
from copyproof.windows import Windows  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--run", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--truth-label", default="truth")
ap.add_argument("--window", type=int, default=1000)
a = ap.parse_args()

run = Path(a.run)
sc = Scan(str(run / "scan"))
ti = sc.reads.index(a.truth_label)
w = Windows.load(str(run / f"windows.w{a.window}.npz"))
model = Model.from_json((run / "model.json").read_text())
summary = json.loads((run / "summary.json").read_text())
res = {"run": str(run), "assembly_bases": int(sc.contigs["length"].sum()), "sequences": int(len(sc.contigs)),
       "depth_per_copy": summary["depth_per_copy"], "assessable_fraction": summary["assessable_fraction"]}

# ---- truth labels: independent of the calls, so computed once per scan and kept ---------------------
static_path = Path(f"{a.out}.truth_static.json")
raw_path = Path(f"{a.out}.truth_cn_segments_raw.tsv")
if static_path.is_file() and raw_path.is_file():
    static = json.loads(static_path.read_text())
    tseg = pd.read_csv(raw_path, sep="\t")
else:
    tw = truth_windows(sc, ti, W=w.W, cap=w.cap)
    tseg = truth_segments(tw, sc.contigs)
    tseg.to_csv(raw_path, sep="\t", index=False)
    del tw
    # read k-mers missing from the assembly: true if the truth genome holds them
    from copyproof.missing import MissingParams, solid_missing  # noqa: E402

    kmers, counts, thr, cand = solid_missing(sc, model.lam, MissingParams())
    in_truth = counts[:, ti] > 0
    truth_absent = int(sc.meta["reads"][ti]["absent"])  # distinct tracked truth k-mers that the assembly lacks
    static = {"missing_kmers": {
        "threshold": thr,
        "solid_missing_tracked": int(len(kmers)),
        "solid_missing_in_truth": int(in_truth.sum()),
        "precision_all_solid": round(float(in_truth.mean()), 4) if len(kmers) else None,
        "candidates_tracked": int(cand.sum()),
        "candidates_in_truth": int((cand & in_truth).sum()),
        "precision_candidates": round(float((cand & in_truth).sum() / max(cand.sum(), 1)), 4),
        "truth_kmers_absent_from_assembly_tracked": truth_absent,
        "recall_all_solid": round(float(in_truth.sum() / max(truth_absent, 1)), 4),
        "recall_candidates": round(float((cand & in_truth).sum() / max(truth_absent, 1)), 4),
        "estimated_missing_truth_kmers": truth_absent * sc.sample,
    }}
    # assembly k-mers without read support: a true error if the truth lacks the k-mer
    tp = fp = fn = tot = 0
    for ci in range(len(sc.contigs)):
        pos, m, c = sc.records(ci)
        mt = m.sum(axis=1)
        el = mt <= w.cap
        no_reads = c[el, 0] == 0
        no_truth = c[el, ti] == 0
        tp += int((no_reads & no_truth).sum())
        fp += int((no_reads & ~no_truth).sum())
        fn += int((~no_reads & no_truth).sum())
        tot += int(el.sum())
    static["unsupported_kmers"] = {
        "positions": tot, "wrong_in_truth": tp + fn, "unsupported_by_reads": tp + fp, "both": tp,
        "precision": round(tp / max(tp + fp, 1), 4), "recall": round(tp / max(tp + fn, 1), 4),
        "wrong_kmer_rate": (tp + fn) / max(tot, 1),
    }
    static_path.write_text(json.dumps(static, indent=1))

# ---- copy number -------------------------------------------------------------------------------
mask = pd.read_csv(run / "cannot_assess.bed", sep="\t", names=["contig", "start", "end"]) if (run / "cannot_assess.bed").stat().st_size else pd.DataFrame(columns=["contig", "start", "end"])
mby = {c: g[["start", "end"]].values for c, g in mask.groupby("contig")}


def masked_frac(c, s, e):
    m = mby.get(c)
    if m is None:
        return 0.0
    ov = np.minimum(m[:, 1], e) - np.maximum(m[:, 0], s)
    return float(ov[ov > 0].sum()) / max(e - s, 1)


tseg["masked"] = [masked_frac(c, s, e) for c, s, e in zip(tseg["contig"], tseg["start"], tseg["end"])]
tseg["len"] = tseg["end"] - tseg["start"]
segs = pd.read_csv(run / "cn_segments.tsv", sep="\t")
# What a user sees: the copy-number issues, with their confidence. A segment that is not an issue
# (a depth anomaly, or ambiguous without support) is scored as "not reported".
iss = pd.read_csv(run / "issues.tsv", sep="\t")
iss = iss[(iss["channel"] == "S1") & iss["class"].isin(["COLLAPSE", "EXPANSION", "UNSUPPORTED"])].copy()
iss["direction"] = np.where(iss["class"] == "COLLAPSE", "gain", "loss")
iss["len"] = iss["end"] - iss["start"]
iss, tseg = match_segments(iss, tseg)
rest = segs[~segs["kind"].isin(["CN", "ABSENT"])].copy()
key = set(zip(iss["contig"], iss["start"], iss["end"]))
rest = rest[[k not in key for k in zip(rest["contig"], rest["start"], rest["end"])]]
rest["len"] = rest["end"] - rest["start"]
rest, _ = match_segments(rest, tseg) if len(rest) else (rest.assign(true=False), None)
tseg.to_csv(f"{a.out}.truth_cn_segments.tsv", sep="\t", index=False)
iss.to_csv(f"{a.out}.cn_issues_scored.tsv", sep="\t", index=False)
high = iss[iss["confidence"] == "high"]
_, t_high = match_segments(high, tseg.drop(columns=["call_overlap", "found"])) if len(high) else (None, tseg.assign(found=False))
tseg["found_high"] = t_high["found"].values

bins = [0, 5_000, 20_000, 100_000, 10**12]
names = ["2-5 kb", "5-20 kb", "20-100 kb", ">100 kb"]
cn = {"truth_segments": int(len(tseg)), "truth_bp": int(tseg["len"].sum()), "issues": int(len(iss)), "issues_bp": int(iss["len"].sum())}
rows = []
for direction in ("gain", "loss"):
    t = tseg[tseg["direction"] == direction]
    t_ok = t[t["masked"] < 0.5]
    c = iss[iss["direction"] == direction]
    for conf in ("high", "low"):
        ck = c[c["confidence"] == conf]
        cn[f"{direction}_{conf}_issues"] = int(len(ck))
        cn[f"{direction}_{conf}_true"] = int(ck["true"].sum())
    for kind in ("CN_OR_DEPTH", "DEPTH", "IRREGULAR"):
        rk = rest[(rest["direction"] == direction) & (rest["kind"] == kind)]
        cn[f"{direction}_{kind}_not_reported"] = int(len(rk))
        cn[f"{direction}_{kind}_not_reported_true"] = int(rk["true"].sum())
    cn[f"{direction}_truth"] = int(len(t))
    cn[f"{direction}_truth_assessable"] = int(len(t_ok))
    cn[f"{direction}_found_assessable"] = int(t_ok["found"].sum())
    cn[f"{direction}_found_high_assessable"] = int(t_ok["found_high"].sum())
    cn[f"{direction}_truth_bp_assessable"] = int(t_ok["len"].sum())
    cn[f"{direction}_found_bp_assessable"] = int(t_ok.loc[t_ok["found"], "len"].sum())
    for lo, hi, nm in zip(bins[:-1], bins[1:], names):
        tt = t_ok[(t_ok["len"] >= lo) & (t_ok["len"] < hi)]
        cc = c[(c["len"] >= lo) & (c["len"] < hi)]
        ch = cc[cc["confidence"] == "high"]
        rows.append({"direction": direction, "size": nm, "truth_assessable": len(tt), "found": int(tt["found"].sum()),
                     "recall": round(float(tt["found"].mean()), 3) if len(tt) else np.nan,
                     "found_high": int(tt["found_high"].sum()), "issues": len(cc), "issues_true": int(cc["true"].sum()),
                     "precision": round(float(cc["true"].mean()), 3) if len(cc) else np.nan, "high": len(ch),
                     "high_true": int(ch["true"].sum()),
                     "precision_high": round(float(ch["true"].mean()), 3) if len(ch) else np.nan})
res["copy_number"] = cn
by_size = pd.DataFrame(rows)
by_size.to_csv(f"{a.out}.cn_by_size.tsv", sep="\t", index=False)
res.update(static)
# read sequence the assembly lacks, localised (only when the run had the reads): a cluster is true
# if the truth genome holds most of its missing k-mers
mc = run / "missing_clusters.tsv"
if mc.is_file():
    c = pd.read_csv(mc, sep="\t")
    col = f"frac_in_{a.truth_label}"
    if col in c.columns and len(c):
        c["true"] = c[col] >= 0.5
        c["strong"] = (c["n_reads"] >= 6) & (c["tracked_kmers"] >= 5)
        res["missing_clusters"] = {
            f"{kind}|{'high' if strong else 'low'}": {"n": int(len(g)), "true": int(g["true"].sum())}
            for (kind, strong), g in c.groupby(["kind", "strong"])
        }
res["issues_by_class"] = {f"{c}|{k}": int(n) for (c, k), n in pd.read_csv(run / "issues.tsv", sep="\t").groupby(["class", "confidence"]).size().items()}
res["kmer_support"] = summary.get("kmer_support", {})
Path(f"{a.out}.truth_eval.json").write_text(json.dumps(res, indent=1))
pd.set_option("display.width", 200)
print(json.dumps(res, indent=1))
print(by_size.to_string(index=False))
