#!/usr/bin/env python3
"""The tables of this record, from the run directories of the held-out test.

  tables.py DATA OUT_DIR

DATA is the data root (the runs are under DATA/work/rung6). Each part is written when its inputs
exist, so the script is run after the calibration (it then writes predictions.tsv, before NA19240
is run), after the scored run, and again at the end. Nothing here decides anything: the scores are
computed by scripts/eval_assemblies.py and scripts/eval_planted.py; this script copies their
tables, sets them beside the predictions and adds the descriptive tables.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

here = Path(__file__).resolve()
sys.path.insert(0, str(here.parents[2] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402

D, out = Path(sys.argv[1]), Path(sys.argv[2])
W = D / "work" / "rung6"
CAL = W / "hg002_calibration" / "hprc_year2"
A = D / "hprc_r2" / "NA19240"
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 300)


def copy(src, name) -> bool:
    src = Path(src)
    if src.is_file():
        shutil.copyfile(src, out / name)
        return True
    return False


def read_bed(path) -> pd.DataFrame:
    """A published QC BED (contig, start, end, label, ...), with or without a track line."""
    rows = [ln.rstrip("\n").split("\t")[:4] for ln in open(path)
            if ln.strip() and not ln.startswith(("#", "track", "browser"))]
    return pd.DataFrame(rows, columns=["contig", "start", "end", "label"]).astype({"start": int, "end": int})


def table(path) -> pd.DataFrame:
    """A table of the run; the hapmer audit writes an empty file when it has nothing to list."""
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def pick(scores: pd.DataFrame, labels: str, direction: str, size: str):
    r = scores[(scores["labels"] == labels) & (scores["direction"] == direction) & (scores["size"] == size)]
    return r.iloc[0]


def planted_sums(t: pd.DataFrame) -> dict:
    """The four sums of entry p from a table of scripts/eval_planted.py."""
    def s(kinds, signal, lo):
        x = t[t["kind"].isin(kinds) & (t["signal"] == signal) & (t["size"] >= lo)]
        return int(x["found"].sum()), int(x["planted"].sum())
    return {"3a": s(["dup", "extra"], "EXPANSION", 5000), "3b": s(["del"], "COLLAPSE", 5000),
            "3c": s(["del", "swap"], "MISSING_HAPLOTYPE", 5000), "3d": s(["swap"], "PHASE_SWITCH", 0)}


WHAT = {
    "1": "high-confidence copy-number issues confirmed by both other assemblies",
    "2 gain": "assessable discordances of 5 kb and more found: the other assemblies hold more (collapse)",
    "2 loss": "assessable discordances of 5 kb and more found: the other assemblies hold fewer (surplus copy)",
    "3a": "planted false duplications found, from 5 kb",
    "3b": "planted deletions found as a collapse on the homologue, from 5 kb",
    "3c": "planted deleted or overwritten stretches found as a missing haplotype, from 5 kb",
    "3d": "planted swaps found as a phase switch",
}

# ---- A. the calibration on HG002 and the predictions ------------------------------------------
if copy(f"{CAL}.eval.cn_scores.tsv", "hg002_cn_scores.tsv"):
    copy(f"{CAL}.eval.calibration.tsv", "hg002_proxy_against_truth.tsv")
    copy(f"{CAL}.eval.missing_clusters.tsv", "hg002_missing_clusters.tsv")
    copy(f"{CAL}.eval.unread_kmers.tsv", "hg002_unread_kmers.tsv")
    copy(f"{CAL}.eval_pair2.cn_scores.tsv", "hg002_pair2_cn_scores.tsv")
    copy(f"{CAL}.eval_pair2.calibration.tsv", "hg002_pair2_proxy_against_truth.tsv")
    copy(D / "work" / "plant" / "run35_trio.eval.sensitivity.tsv", "hg002_planted_with_parents.sensitivity.tsv")
    cs = pd.read_csv(out / "hg002_cn_scores.tsv", sep="\t")
    rows = []
    r = pick(cs, "both", "both directions", "all")
    rows.append(("1", int(r["high_confirmed"]), int(r["high"]), "hg002_cn_scores.tsv"))
    for d in ("gain", "loss"):
        r = pick(cs, "both", d, "5 kb and more")
        rows.append((f"2 {d}", int(r["found"]), int(r["assessable"]), "hg002_cn_scores.tsv"))
    ps = planted_sums(pd.read_csv(out / "hg002_planted_with_parents.sensitivity.tsv", sep="\t"))
    rows += [(k, v[0], v[1], "hg002_planted_with_parents.sensitivity.tsv") for k, v in ps.items()]
    pred = pd.DataFrame(rows, columns=["score", "predicted_count", "predicted_of", "source"])
    pred.insert(1, "what", pred["score"].map(WHAT))
    pred.to_csv(out / "predictions.tsv", sep="\t", index=False)
    print("predictions (HG002):")
    print(pred.to_string(index=False))
    # the calibration run next to the truth's verdict on the same issues
    cal = pd.read_csv(out / "hg002_proxy_against_truth.tsv", sep="\t")
    print(cal[cal["what"] == "issues"].to_string(index=False))

# ---- B. NA19240: the scored run and the run with long reads -----------------------------------
runs = {name: W / f"NA19240_{name}" for name in ("short", "long")}
have = {name: Path(f"{p}.eval.cn_scores.tsv").is_file() for name, p in runs.items()}
for name, p in runs.items():
    if not have[name]:
        continue
    copy(f"{p}.eval.cn_scores.tsv", f"na19240_{name}_cn_scores.tsv")
    copy(f"{p}.eval.missing_clusters.tsv", f"na19240_{name}_missing_clusters.tsv")
    copy(f"{p}.eval.unread_kmers.tsv", f"na19240_{name}_unread_kmers.tsv")
    copy(p / "summary.json", f"summary_{name}.json")
    copy(p / "run.json", f"run_{name}.json")

if have["short"]:
    # issues by class, beside the run with long reads and beside HG01943 (rung 1)
    tabs = []
    for name, p in runs.items():
        if (p / "issues.tsv").is_file():
            iss = pd.read_csv(p / "issues.tsv", sep="\t")
            tabs.append(iss.groupby(["class", "confidence", "role"]).size().rename(f"NA19240_{name}"))
    old = out.parent / "2026-11-13_hg01943_trio" / "issues_by_class.tsv"
    if old.is_file():
        tabs.append(pd.read_csv(old, sep="\t").set_index(["class", "confidence", "role"])["n"].rename("HG01943_long"))
    pd.concat(tabs, axis=1).fillna(0).astype(int).sort_index().reset_index().to_csv(out / "issues_by_class.tsv", sep="\t", index=False)

    # the high-confidence copy-number issues that the two other assemblies do not confirm together
    p = runs["short"]
    ci = pd.read_csv(f"{p}.eval.cn_issues.tsv", sep="\t")
    iss = pd.read_csv(p / "issues.tsv", sep="\t").set_index("id")
    un = ci[(ci["confidence"] == "high") & ~ci["confirmed_both"]].copy()
    ev = iss.loc[un["id"], "evidence"]
    for key in ("delta", "copies_by_multiplicity", "windows"):
        un[key] = ev.str.extract(rf"{key}=([^;]+)")[0].values
    beds = []
    for f in sorted((A / "qc").glob("*.bed")):
        b = read_bed(f)
        b["tool"] = "flagger" if "flagger" in f.name else "nucflag"
        beds.append(b)
    beds = pd.concat(beds, ignore_index=True)
    qc = {}
    for which in ("flagger", "nucflag"):
        b = beds[(beds["tool"] == which) & (beds["label"] != "NNN")]
        qc[which] = IntervalSet.from_rows(zip(b["contig"], b["start"], b["end"]))
    for which, s in qc.items():
        un[f"at_{which}_call"] = [s.overlap(c, a, b) > 0 for c, a, b in zip(un["contig"], un["start"], un["end"])]
    refs = json.loads(Path(f"{p}.eval.json").read_text())["refs"]
    un.sort_values([f"overlap_{r}" for r in refs], ascending=False).to_csv(out / "unconfirmed_high.tsv", sep="\t", index=False)

    # overlap with HPRC's published QC, both ways (the run with long reads when there is one)
    q = runs["long"] if (runs["long"] / "issues.tsv").is_file() else p
    iq = pd.read_csv(q / "issues.tsv", sep="\t")
    for which, s in qc.items():
        iq[which] = [s.overlap(c, a, b) > 0 for c, a, b in zip(iq["contig"], iq["start"], iq["end"])]
    g = iq.groupby(["class", "confidence"])
    t = g.agg(issues=("id", "size"), at_flagger_call=("flagger", "sum"), at_nucflag_call=("nucflag", "sum")).reset_index()
    t.insert(0, "run", q.name)
    t.to_csv(out / "issues_at_published_calls.tsv", sep="\t", index=False)
    everything = IntervalSet.from_rows(zip(iq["contig"], iq["start"], iq["end"]))
    cn = iq[iq["class"].isin(["COLLAPSE", "EXPANSION", "UNSUPPORTED"])]
    cn_set = IntervalSet.from_rows(zip(cn["contig"], cn["start"], cn["end"]))
    b = beds.copy()
    b["tool"] = b["tool"].replace({"flagger": "flagger_hifi_conservative"})
    b["any_issue"] = [everything.overlap(c, a, e) > 0 for c, a, e in zip(b["contig"], b["start"], b["end"])]
    b["cn_issue"] = [cn_set.overlap(c, a, e) > 0 for c, a, e in zip(b["contig"], b["start"], b["end"])]
    b["bp"] = b["end"] - b["start"]
    b.groupby(["tool", "label"]).agg(calls=("bp", "size"), bp=("bp", "sum"), at_any_issue=("any_issue", "sum"),
                                     at_copy_number_issue=("cn_issue", "sum")).reset_index().to_csv(
        out / "published_calls_at_issues.tsv", sep="\t", index=False)

    # hapmer channel and, with long reads, the joins and the long-read verdicts
    rows = []
    for name, r in runs.items():
        sp = r / "summary.json"
        if not sp.is_file():
            continue
        s = json.loads(sp.read_text())
        row = {"run": r.name, "depth_per_copy": s["depth_per_copy"], "assessable_fraction": s["assessable_fraction"],
               "issues": s["issues"], "read_kmer_qv": s["kmer_support"].get("qv"), **{f"hapmers_{k}": v for k, v in s["hapmers"].items()}}
        h = table(r / "hapmer_switches.tsv") if (r / "hapmer_switches.tsv").is_file() else pd.DataFrame()
        if len(h):
            row.update(switch_segments_10_sites=int((h["sites_segment_parent"] >= 10).sum()), switch_bp=int((h["end"] - h["start"]).sum()))
        h = table(r / "hapmer_pieces.tsv") if (r / "hapmer_pieces.tsv").is_file() else pd.DataFrame()
        row["pieces"] = len(h)
        if len(h):
            row.update({f"pieces_{k}": int(v) for k, v in h["call"].value_counts().items()})
        j = s.get("joins", {})
        row.update({f"joins_{k}": v for k, v in j.items() if not isinstance(v, dict)})
        for k in ("unsupported_runs", "copy_number_losses"):
            row.update({f"{k}_{kk}": vv for kk, vv in j.get(k, {}).items()})
        row.update({f"seconds_{k[:-2]}": v for k, v in s["timings"].items()})
        rows.append(row)
    # each value as the run wrote it (a count stays a count), nothing where a run has no such entry
    t = pd.DataFrame([{k: (v if k == "run" else "" if v is None else str(v)) for k, v in row.items()} for row in rows])
    t = t.set_index("run").T.fillna("")
    t.reset_index().rename(columns={"index": "what"}).to_csv(out / "runs.tsv", sep="\t", index=False)

# ---- C. planted errors ------------------------------------------------------------------------
P = W / "plant"
if copy(f"{P}/run.eval.sensitivity.tsv", "planted_sensitivity.tsv"):
    copy(P / "plant1.edits.tsv", "planted_edits.tsv")
    copy(f"{P}/run.eval.new_elsewhere.tsv", "planted_new_elsewhere.tsv")
    copy(f"{P}/run.eval.breakpoints.tsv", "planted_breakpoints.tsv")
    copy(P / "run" / "summary.json", "summary_planted.json")

# ---- the added run: the hapmer audit on the Verkko assembly, held against HPRC's haplotypes ----
V = W / "NA19240_verkko"
if copy(f"{V}.hapmers.summary.tsv", "verkko_hapmers_summary.tsv"):
    copy(f"{V}.hapmers.calls.tsv", "verkko_hapmers_calls.tsv")
    copy(f"{V}.hapmers.sequences.tsv", "verkko_hapmers_sequences.tsv")
    copy(f"{V}.hapmers.other_parent_bins.tsv", "verkko_hapmers_other_parent_bins.tsv")
    copy(V / "summary.json", "summary_verkko.json")
    copy(V / "run.json", "run_verkko.json")

# ---- what was counted: the reads of the trio and the assemblies used as label sets ----
rows = []
for name, what in (("NA19240.nygc", "reads of NA19240 (NYGC 30x CRAM)"), ("NA19238.nygc", "reads of the mother, NA19238"),
                   ("NA19239.nygc", "reads of the father, NA19239"),
                   ("NA19240.hgsvc3_verkko.exact", "HGSVC3 Verkko assembly, every file"),
                   ("NA19240.hgsvc3_hifiasm.exact", "HGSVC3 hifiasm assembly, both haplotypes"),
                   ("NA19240.hprc_r2_pat.exact", "HPRC release 2, paternal haplotype"),
                   ("NA19240.hprc_r2_mat.exact", "HPRC release 2, maternal haplotype")):
    f = D / "kmerdb" / f"{name}.json"
    if not f.is_file():
        continue
    j = json.loads(f.read_text())
    h = j["header"]
    rows.append({"database": f"{name}.cpkdb", "what": what, "every_sighting_counted": h["exact"], "records": h["records"],
                 "bases": h["bases"], "kmers": h["kmers"], "tracked_kmers": h["tracked"], "distinct_tracked_kmers_kept": h["n_entries"],
                 "kept_from_count": h["min_count"], "reads_cut_at_adapter": h["adapter_trimmed_reads"],
                 "seconds": j["seconds"], "peak_rss_gb": j["peak_rss_gb"]})
if rows:
    pd.DataFrame(rows).to_csv(out / "kmer_databases.tsv", sep="\t", index=False)

# ---- score 1 as counts: confirmed and not confirmed, for every label set (entry p, second addition) ----
rows = []
for run, name in (("HG002 hprc_year2, judged by v2_trio and hprc_year1", "hg002_cn_scores.tsv"),
                  ("HG002 hprc_year2, judged by verkko_asm_thic and hprc_year2_v2_thic", "hg002_pair2_cn_scores.tsv"),
                  ("NA19240, short reads of the trio", "na19240_short_cn_scores.tsv"),
                  ("NA19240, with long reads", "na19240_long_cn_scores.tsv")):
    if not (out / name).is_file():
        continue
    t = pd.read_csv(out / name, sep="\t")
    t = t[t["size"] == "all"]
    for r in t.itertuples(index=False):
        rows.append({"run": run, "judged_by": r.labels, "direction": r.direction, "high_confidence_issues": int(r.high),
                     "confirmed": int(r.high_confirmed), "not_confirmed": int(r.high - r.high_confirmed),
                     "share_confirmed": round(r.high_confirmed / r.high, 3) if r.high else np.nan})
if rows:
    pd.DataFrame(rows).to_csv(out / "score1_counts.tsv", sep="\t", index=False)

# ---- D. the verdicts: each score beside its prediction ----------------------------------------
if (out / "predictions.tsv").is_file() and have["short"]:
    pred = pd.read_csv(out / "predictions.tsv", sep="\t", dtype={"score": str}).set_index("score")
    cs = pd.read_csv(out / "na19240_short_cn_scores.tsv", sep="\t")
    got = {}
    r = pick(cs, "both", "both directions", "all")
    got["1"] = (int(r["high_confirmed"]), int(r["high"]))
    for d in ("gain", "loss"):
        r = pick(cs, "both", d, "5 kb and more")
        got[f"2 {d}"] = (int(r["found"]), int(r["assessable"]))
    if (out / "planted_sensitivity.tsv").is_file():
        got.update(planted_sums(pd.read_csv(out / "planted_sensitivity.tsv", sep="\t")))
    rows = []
    for score, (k, n) in got.items():
        k0, n0 = int(pred.at[score, "predicted_count"]), int(pred.at[score, "predicted_of"])
        p = float(fisher_exact([[k, n - k], [k0, n0 - k0]], alternative="less")[1]) if n and n0 else float("nan")
        rows.append({"score": score, "what": WHAT[score], "predicted": f"{k0} of {n0}", "predicted_share": round(k0 / n0, 3) if n0 else np.nan,
                     "NA19240": f"{k} of {n}", "NA19240_share": round(k / n, 3) if n else np.nan, "p_lower_one_sided": round(p, 4),
                     "holds": "no" if p < 0.05 else ("yes" if n else "no count")})
    v = pd.DataFrame(rows)
    v.to_csv(out / "verdicts.tsv", sep="\t", index=False)
    print(v.to_string(index=False))
