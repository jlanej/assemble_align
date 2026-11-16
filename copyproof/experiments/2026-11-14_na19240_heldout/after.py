#!/usr/bin/env python3
"""Tables added after the scored run of NA19240 had been seen. They decide nothing.

`tables.py` writes the scores and the verdicts and was fixed before the run. This script was
written afterwards, to describe what the scores leave open: which issues the two other assemblies
do not confirm and how, which discordances no issue covers and what the copy-number channel holds
there, whether the missing-allele calls that neither assembly holds are in the parents' reads, and
how the overlap with HPRC's published calls compares with the overlap of issues placed at random.

  after.py DATA OUT_DIR

It also rewrites the data root as DATA in the copies of the runs' parameters (run_*.json).
Writes, each when its inputs exist:
  after_confirmed_by.tsv            high-confidence copy-number issues by class and by which assembly confirms
  after_missed_discordances.tsv     assessable discordances by `both` of 5 kb and more that no issue covers
  after_missed_summary.tsv          the same, counted by what the copy-number channel holds there
  after_reads_as_referee.tsv        every label set's discordances of 5 kb and more: is the read depth nearer to the release or to the others?
  after_missing_alleles.tsv         strong missing-allele calls: held by the other assemblies, in the parents' reads
  after_missing_alleles_by_parent.tsv   the same, summed by the parents' reads
  after_issues_at_calls.tsv         issues at HPRC's calls, beside the same issues placed at random on their sequences
  after_issues_by_call_label.tsv    the same issues by the label of the call they lie at
  after_published_calls.tsv         HPRC's calls inside and outside the run's mask, and at issues
  after_long_reads_at_copy_number_issues.tsv, after_long_reads_summary.tsv, after_long_reads_score1.tsv
                                    every copy-number issue of the short run followed into the run with long reads
  after_against_hg002.tsv           NA19240 beside both pairs of HG002 assemblies, and beside each judge alone
  after_planted_totals.tsv          issues on the planted genome at edits and elsewhere
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

here = Path(__file__).resolve()
sys.path.insert(0, str(here.parents[2] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402

D, out = Path(sys.argv[1]), Path(sys.argv[2])
W = D / "work" / "rung6"
RUN = W / "NA19240_short"
CAL = W / "hg002_calibration" / "hprc_year2"
A = D / "hprc_r2" / "NA19240"
CN = ["COLLAPSE", "EXPANSION", "UNSUPPORTED"]
DRAWS, SEED = 200, 1


def read_bed(path) -> pd.DataFrame:
    rows = [ln.rstrip("\n").split("\t")[:4] for ln in open(path)
            if ln.strip() and not ln.startswith(("#", "track", "browser"))]
    return pd.DataFrame(rows, columns=["contig", "start", "end", "label"]).astype({"start": int, "end": int})


def iset(df) -> IntervalSet:
    return IntervalSet.from_rows(zip(df["contig"], df["start"], df["end"]))


def hits(s: IntervalSet, df) -> np.ndarray:
    return np.array([s.overlap(c, a, b) > 0 for c, a, b in zip(df["contig"], df["start"], df["end"])], dtype=bool)


if not Path(f"{RUN}.eval.cn_issues.tsv").is_file():
    sys.exit("the scored run has not been scored yet")

# The copies of each run's parameters name the inputs by their place on the machine that ran them.
# The documents name them under DATA; so do the copies from here on.
for f in out.glob("run_*.json"):
    text = f.read_text()
    if str(D) in text:
        f.write_text(text.replace(str(D), "DATA"))

beds = []
for f in sorted((A / "qc").glob("*.bed")):
    b = read_bed(f)
    b["tool"] = "flagger" if "flagger" in f.name else "nucflag"
    beds.append(b)
beds = pd.concat(beds, ignore_index=True)
qc = {t: iset(beds[(beds["tool"] == t) & (beds["label"] != "NNN")]) for t in ("flagger", "nucflag")}
mask = IntervalSet.from_bed(str(RUN / "cannot_assess.bed"))
issues = pd.read_csv(RUN / "issues.tsv", sep="\t")
refs = json.loads(Path(f"{RUN}.eval.json").read_text())["refs"]
ra, rb = refs

# ---- 1. the high-confidence copy-number issues, by which assembly confirms ---------------------
ci = pd.read_csv(f"{RUN}.eval.cn_issues.tsv", sep="\t")
hi = ci[ci["confidence"] == "high"].copy()
for t, s in qc.items():
    hi[t] = hits(s, hi)
ca, cb, both = hi[f"confirmed_{ra}"], hi[f"confirmed_{rb}"], hi["confirmed_both"]
hi["confirmed_by"] = np.select(
    [both, ca & cb, ca, cb], ["both, with the same number of copies", f"{ra} and {rb}, with different numbers of copies",
                              f"{ra} only", f"{rb} only"], "neither")
t1 = hi.groupby(["class", "direction", "confirmed_by"]).agg(
    issues=("id", "size"), bp=("len", "sum"), at_flagger_call=("flagger", "sum"), at_nucflag_call=("nucflag", "sum")).reset_index()
tot = hi.groupby("confirmed_by").agg(issues=("id", "size"), bp=("len", "sum"), at_flagger_call=("flagger", "sum"),
                                     at_nucflag_call=("nucflag", "sum")).reset_index()
tot.insert(0, "direction", "both directions")
tot.insert(0, "class", "all")
pd.concat([t1, tot], ignore_index=True).to_csv(out / "after_confirmed_by.tsv", sep="\t", index=False)

# ---- 2. the discordances no issue covers, and what the copy-number channel holds there ---------
seg = pd.read_csv(f"{RUN}.eval.segments.tsv", sep="\t")
d = seg[(seg["labels"] == "both") & (seg["len"] >= 5000) & (seg["masked"] < 0.5)].copy()
cs = pd.read_csv(RUN / "cn_segments.tsv", sep="\t")
kinds = sorted(cs["kind"].unique())
by_kind = {k: iset(cs[cs["kind"] == k]) for k in kinds}
cn_issue_any = iset(issues[issues["class"].isin(CN)])
for k in kinds:
    d[f"bp_{k}"] = [by_kind[k].overlap(c, a, b) for c, a, b in zip(d["contig"], d["start"], d["end"])]
d["bp_any_cn_issue"] = [cn_issue_any.overlap(c, a, b) for c, a, b in zip(d["contig"], d["start"], d["end"])]
for t, s in qc.items():
    d[f"at_{t}_call"] = hits(s, d)
for r in refs:
    own = iset(seg[(seg["labels"] == r)])
    d[f"bp_discordant_by_{r}"] = [own.overlap(c, a, b) for c, a, b in zip(d["contig"], d["start"], d["end"])]


def holds(row) -> str:
    """The segment kind that covers most of the discordance, or none when under a tenth is covered."""
    best = max(kinds, key=lambda k: row[f"bp_{k}"])
    return best if row[f"bp_{best}"] >= 0.1 * row["len"] else "no segment"


d["copy_number_channel_holds"] = d.apply(holds, axis=1)

# read depth over the discordance, relative to what the release's own copy numbers predict
# (1 when the reads agree with the release; no GC correction, so a few percent either way mean nothing)
z = np.load(RUN / "windows.w1000.npz")
n_arr, S_arr, width = z["n"], z["S"], int(z["W"])
win0 = pd.read_csv(RUN / "windows.w1000.contigs.tsv", sep="\t").set_index("name")["win_start"]
depth = json.loads((RUN / "summary.json").read_text())["depth_per_copy"]
mult = np.arange(n_arr.shape[1])


def reads_over_release(contig, start, end) -> float:
    i, j = int(win0[contig]) + start // width, int(win0[contig]) + -(-end // width)
    expected = depth * float((n_arr[i:j] * mult).sum())
    return round(float(S_arr[i:j].sum()) / expected, 3) if expected > 0 else float("nan")


d["reads_over_release"] = [reads_over_release(c, a, b) for c, a, b in zip(d["contig"], d["start"], d["end"])]
d["others_over_release"] = [round(1 + md / max(float((n_arr[int(win0[c]) + a // width:int(win0[c]) + -(-b // width)] * mult).sum())
                                                 / max(float(n_arr[int(win0[c]) + a // width:int(win0[c]) + -(-b // width)].sum()), 1), 1e-9), 3)
                            for c, a, b, md in zip(d["contig"], d["start"], d["end"], d["mean_delta"])]
# the reads as referee, for every label set: over each assessable discordance of 5 kb and more, is the
# read depth nearer to what the release holds (1) or to what the other assemblies hold?
ref = seg[(seg["len"] >= 5000) & (seg["masked"] < 0.5)].copy()
ref["reads_over_release"] = [reads_over_release(c, a, b) for c, a, b in zip(ref["contig"], ref["start"], ref["end"])]
ref["others_over_release"] = [round(1 + md / max(float((n_arr[int(win0[c]) + a // width:int(win0[c]) + -(-b // width)] * mult).sum())
                                                   / max(float(n_arr[int(win0[c]) + a // width:int(win0[c]) + -(-b // width)].sum()), 1), 1e-9), 3)
                              for c, a, b, md in zip(ref["contig"], ref["start"], ref["end"], ref["mean_delta"])]
ref["reads_nearer_to"] = np.where((ref["reads_over_release"] - 1).abs() <= (ref["reads_over_release"] - ref["others_over_release"]).abs(),
                                  "the release", "the other assemblies")
ref.groupby(["labels", "direction", "found", "reads_nearer_to"]).agg(
    discordances=("len", "size"), bp=("len", "sum"), median_reads_over_release=("reads_over_release", "median"),
    median_others_over_release=("others_over_release", "median")).reset_index().to_csv(out / "after_reads_as_referee.tsv", sep="\t", index=False)

chain_path = W / "plant" / "hap_to_hap.chain"
if chain_path.is_file():
    sys.path.insert(0, str(here.parents[2] / "scripts"))
    from plant_errors import Chain  # noqa: E402

    chain = Chain(str(chain_path))
    lifted = [chain.lift_interval(c, a, b) for c, a, b in zip(d["contig"], d["start"], d["end"])]
    d["homologue"] = [f"{x[0]}:{x[1]}-{x[2]}" if x else "" for x in lifted]
    same = {dr: iset(d[d["direction"] == dr]) for dr in ("gain", "loss")}
    notfound = {dr: iset(d[(d["direction"] == dr) & ~d["found"]]) for dr in ("gain", "loss")}
    d["homologue_is_discordant_too"] = [bool(x) and same[dr].overlap(*x) > 0 for x, dr in zip(lifted, d["direction"])]
    d["homologue_is_not_found_either"] = [bool(x) and notfound[dr].overlap(*x) > 0 for x, dr in zip(lifted, d["direction"])]
missed = d[~d["found"]].sort_values(["direction", "len"], ascending=[True, False])
missed.to_csv(out / "after_missed_discordances.tsv", sep="\t", index=False)
t2 = d.groupby(["direction", "found", "copy_number_channel_holds"]).agg(
    discordances=("len", "size"), bp=("len", "sum"), median_reads_over_release=("reads_over_release", "median"),
    median_others_over_release=("others_over_release", "median"), at_flagger_call=("at_flagger_call", "sum"),
    at_nucflag_call=("at_nucflag_call", "sum")).reset_index()
t2.to_csv(out / "after_missed_summary.tsv", sep="\t", index=False)

# ---- 3. strong missing-allele calls: the other assemblies and the parents' reads ---------------
rows, by_parent = [], []
for name, run, labels, truth in (("NA19240", RUN, refs, None), ("HG002 hprc_year2", CAL, ["v2_trio", "hprc_year1"], "truth")):
    f = run / "missing_clusters.tsv"
    if not f.is_file():
        continue
    c = pd.read_csv(f, sep="\t")
    c = c[(c["kind"] == "MISSING_ALLELE") & (c["n_reads_locus"] >= 6) & (c["tracked_kmers"] >= 5)]  # strong, as in the issue list
    held = np.where((c[f"frac_in_{labels[0]}"] >= 0.5) & (c[f"frac_in_{labels[1]}"] >= 0.5), "both",
                    np.where((c[f"frac_in_{labels[0]}"] >= 0.5) | (c[f"frac_in_{labels[1]}"] >= 0.5), "one", "neither"))
    parent = np.where((c["frac_in_mother"] >= 0.5) | (c["frac_in_father"] >= 0.5), "yes", "no")
    g = pd.DataFrame({"held": held, "parent": parent, "truth": (c[f"frac_in_{truth}"] >= 0.5).values if truth else False,
                      "count": c["median_count"].values})
    for par, x in g.groupby("parent"):
        by_parent.append({"run": name, "depth_per_copy": json.loads((run / "summary.json").read_text())["depth_per_copy"],
                          "in_a_parents_reads": par, "calls": len(x), "held_by_one_or_both": int((x["held"] != "neither").sum()),
                          "share_held": round(float((x["held"] != "neither").mean()), 3), "read_count_median": float(x["count"].median()),
                          "read_count_quartile_1": float(x["count"].quantile(0.25)), "read_count_quartile_3": float(x["count"].quantile(0.75))})
    per_copy = json.loads((run / "summary.json").read_text())["depth_per_copy"]
    for (h, p), x in g.groupby(["held", "parent"]):
        row = {"run": name, "depth_per_copy": per_copy, "held_by_the_other_assemblies": h, "in_a_parents_reads": p, "calls": len(x),
               "read_count_median": float(x["count"].median()), "read_count_quartile_1": float(x["count"].quantile(0.25)),
               "read_count_quartile_3": float(x["count"].quantile(0.75))}
        if truth:
            row["held_by_truth"] = int(x["truth"].sum())
        rows.append(row)
ma = pd.DataFrame(rows)
ma.to_csv(out / "after_missing_alleles.tsv", sep="\t", index=False)
# the same calls by the parents' reads alone: how many does at least one other assembly hold?
pd.DataFrame(by_parent).to_csv(out / "after_missing_alleles_by_parent.tsv", sep="\t", index=False)

# ---- 4. issues at HPRC's calls, beside the same issues placed at random on their sequences ----
length = pd.read_csv(RUN / "scan.contigs.tsv", sep="\t").set_index("name")["length"]
rng = np.random.default_rng(SEED)
rows = []
for (cls, conf), g in issues.groupby(["class", "confidence"]):
    L = (g["end"] - g["start"]).to_numpy()
    n = length.loc[g["contig"]].to_numpy()
    row = {"class": cls, "confidence": conf, "issues": len(g), "median_bp": int(np.median(L))}
    for t, s in qc.items():
        obs = int(hits(s, g).sum())
        exp = []
        for _ in range(DRAWS):
            a = (rng.random(len(g)) * np.maximum(n - L, 1)).astype(np.int64)
            exp.append(sum(s.overlap(c, x, x + l) > 0 for c, x, l in zip(g["contig"], a, L)))
        row.update({f"at_{t}_call": obs, f"at_{t}_call_if_placed_at_random": round(float(np.mean(exp)), 2),
                    f"at_{t}_call_random_max_of_{DRAWS}": int(np.max(exp))})
    rows.append(row)
pd.DataFrame(rows).to_csv(out / "after_issues_at_calls.tsv", sep="\t", index=False)

# which of HPRC's labels the issues lie at (an issue is counted once under every label it touches)
by_label = {(t, lab): iset(g) for (t, lab), g in beds[beds["label"] != "NNN"].groupby(["tool", "label"])}
rows = []
for (cls, conf), g in issues.groupby(["class", "confidence"]):
    row = {"class": cls, "confidence": conf, "issues": len(g)}
    row.update({f"{t}_{lab}": int(hits(s, g).sum()) for (t, lab), s in by_label.items()})
    rows.append(row)
pd.DataFrame(rows).to_csv(out / "after_issues_by_call_label.tsv", sep="\t", index=False)

# ---- 5. HPRC's calls inside and outside the run's mask, and at issues --------------------------
b = beds.copy()
b["bp"] = b["end"] - b["start"]
b["masked_bp"] = [mask.overlap(c, s, e) for c, s, e in zip(b["contig"], b["start"], b["end"])]
b["mostly_masked"] = b["masked_bp"] >= 0.5 * b["bp"].clip(lower=1)
everything = iset(issues)
cn_set = iset(issues[issues["class"].isin(CN)])
b["any_issue"], b["cn_issue"] = hits(everything, b), hits(cn_set, b)
t5 = b.groupby(["tool", "label", "mostly_masked"]).agg(calls=("bp", "size"), bp=("bp", "sum"), at_any_issue=("any_issue", "sum"),
                                                       at_copy_number_issue=("cn_issue", "sum")).reset_index()
t5.to_csv(out / "after_published_calls.tsv", sep="\t", index=False)

# ---- 6. the run with long reads: what the long-read audit did to the copy-number issues ---------
# scripts/eval_assemblies.py scores issues whose channel is exactly S1. In a run with long reads the
# audit marks the issues it looked at S1+S4, so the scores of that run leave them out. Here every
# copy-number issue of the short run is followed into the run with long reads, by its coordinates,
# with the verdict of the other assemblies that the short run's scoring gave it.
names = ["after_confirmed_by", "after_missed_summary", "after_reads_as_referee", "after_missing_alleles", "after_missing_alleles_by_parent", "after_issues_at_calls",
         "after_issues_by_call_label", "after_published_calls"]
LONG = W / "NA19240_long"
if (LONG / "issues.tsv").is_file():
    li = pd.read_csv(LONG / "issues.tsv", sep="\t")
    lc = li[li["class"].isin(CN)].copy()
    lc["long_reads"] = np.where(lc["summary"].str.contains("long reads carry"), "carry its k-mers",
                                np.where(lc["summary"].str.contains("long reads do not carry"), "do not carry them", "not asked"))
    key = ["contig", "start", "end", "class"]
    m = ci.merge(lc[key + ["confidence", "channel", "long_reads"]], on=key, how="outer", suffixes=("", "_long"), indicator=True)
    assert (m["_merge"] != "right_only").all(), "a copy-number issue of the run with long reads is not in the short run"
    m["confidence_long"] = m["confidence_long"].fillna("not reported")
    m["long_reads"] = m["long_reads"].fillna("not reported")
    m = m.drop(columns=["_merge"]).rename(columns={"confidence": "confidence_short", "channel": "channel_long"})
    m.to_csv(out / "after_long_reads_at_copy_number_issues.tsv", sep="\t", index=False)
    t6 = m.groupby(["direction", "confidence_short", "long_reads", "confidence_long"]).agg(
        issues=("id", "size"), **{f"confirmed_{r}": (f"confirmed_{r}", "sum") for r in refs},
        confirmed_both=("confirmed_both", "sum"), confirmed_either=("confirmed_either", "sum")).reset_index()
    t6.to_csv(out / "after_long_reads_summary.tsv", sep="\t", index=False)
    # score 1 over every copy-number issue that is high-confidence in the run with long reads, whatever its channel
    h = m[m["confidence_long"] == "high"]
    rows = [{"run": "NA19240, with long reads, every channel", "judged_by": r, "direction": dname, "high_confidence_issues": len(x),
             "confirmed": int(x[f"confirmed_{r}"].sum()), "not_confirmed": int((~x[f"confirmed_{r}"]).sum())}
            for r in refs + ["both", "either"]
            for dname, x in (("gain", h[h["direction"] == "gain"]), ("loss", h[h["direction"] == "loss"]), ("both directions", h))]
    pd.DataFrame(rows).to_csv(out / "after_long_reads_score1.tsv", sep="\t", index=False)
    names += ["after_long_reads_summary", "after_long_reads_score1"]

# ---- 7. NA19240 beside both pairs of HG002, and beside each judge alone ------------------------
# Entry p declared the prediction from the first pair and said that the record would set NA19240
# against the second pair too. The test is the one of the verdicts (one-sided, is NA19240 lower?);
# only the rows of the declared pair judged by `both` are verdicts, and tables.py writes those.
from scipy.stats import fisher_exact  # noqa: E402

na = pd.read_csv(f"{RUN}.eval.cn_scores.tsv", sep="\t")


def cell(t, labels, direction, size, k, n):
    r = t[(t["labels"] == labels) & (t["direction"] == direction) & (t["size"] == size)].iloc[0]
    return int(r[k]), int(r[n])


rows = []
for pair, f, (ja, jb) in (("declared pair", "hg002_cn_scores.tsv", ("v2_trio", "hprc_year1")),
                          ("second pair", "hg002_pair2_cn_scores.tsv", ("verkko_asm_thic", "hprc_year2_v2_thic"))):
    if not (out / f).is_file():
        continue
    hg = pd.read_csv(out / f, sep="\t")
    for what, direction, size, k, n in (("high-confidence copy-number issues confirmed", "both directions", "all", "high_confirmed", "high"),
                                        ("high-confidence collapse issues confirmed", "gain", "all", "high_confirmed", "high"),
                                        ("high-confidence surplus-copy issues confirmed", "loss", "all", "high_confirmed", "high"),
                                        ("discordances of 5 kb and more found, others hold more", "gain", "5 kb and more", "found", "assessable"),
                                        ("discordances of 5 kb and more found, others hold fewer", "loss", "5 kb and more", "found", "assessable")):
        for judge, lh, ln in (("both", "both", "both"), ("first assembly alone", ja, ra), ("second assembly alone", jb, rb)):
            k0, n0 = cell(hg, lh, direction, size, k, n)
            k1, n1 = cell(na, ln, direction, size, k, n)
            pv = float(fisher_exact([[k1, n1 - k1], [k0, n0 - k0]], alternative="less")[1]) if n0 and n1 else float("nan")
            rows.append({"hg002_judged_by": pair, "judge": judge, "hg002_label": lh, "na19240_label": ln, "what": what, "HG002": f"{k0} of {n0}",
                         "HG002_share": round(k0 / n0, 3) if n0 else np.nan, "NA19240": f"{k1} of {n1}",
                         "NA19240_share": round(k1 / n1, 3) if n1 else np.nan, "p_lower_one_sided": round(pv, 4)})
pd.DataFrame(rows).to_csv(out / "after_against_hg002.tsv", sep="\t", index=False)
names.append("after_against_hg002")

# ---- 8. the planted genome: issues at edits and elsewhere (scripts/eval_planted.py prints these) ----
rows = []
for name, run, pre in (("NA19240", W / "plant" / "run", W / "plant" / "run.eval"),
                       ("HG002 with parents", D / "work" / "plant" / "run35_trio", D / "work" / "plant" / "run35_trio.eval")):
    if not Path(f"{pre}.elsewhere.tsv").is_file():
        continue
    total = len(pd.read_csv(run / "issues.tsv", sep="\t", usecols=["id"]))
    away = pd.read_csv(f"{pre}.elsewhere.tsv", sep="\t", usecols=["id", "in_baseline"])
    listed = len(pd.read_csv(f"{pre}.new_elsewhere.tsv", sep="\t", usecols=["id"]))
    new = int((~away["in_baseline"]).sum())
    rows.append({"planted_genome": name, "issues": total, "at_planted_edits": total - len(away), "elsewhere": len(away),
                 "also_on_the_unedited_genome": int(away["in_baseline"].sum()), "new": new,
                 "new_near_the_homologue_of_an_edit": new - listed, "new_elsewhere": listed})
if rows:
    pd.DataFrame(rows).to_csv(out / "after_planted_totals.tsv", sep="\t", index=False)
    names.append("after_planted_totals")

pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 300)
for name in names:
    print(f"== {name}")
    print(pd.read_csv(out / f"{name}.tsv", sep="\t").to_string(index=False))
