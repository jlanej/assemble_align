#!/usr/bin/env python3
"""The tables of this record: alleles that the reads carry and a release lacks, sorted by whether
a parent's reads carry them.

  tables.py DATA OUT_DIR

DATA is the data root. Only finished runs under DATA/work and the side files of their read
databases are read; nothing is run and nothing under DATA is written.

Definitions, all taken from code that already exists:
  placement   one row of a run's missing_clusters.tsv: reads that agree on one place of one sequence
  event       placements that hold mostly the same reads (column twin_group); an allele whose
              flanks exist on both haplotypes is one event with two placements
  strong      6 reads or more and 5 tracked k-mers or more (src/copyproof/calls.py: what makes a
              missing allele a high-confidence issue)
  in a parent's reads   half or more of the placement's tracked k-mers are in the mother's read
              database, or half or more in the father's (the cut of after.py in the record of the
              held-out test). A read database keeps a k-mer from its second sighting.
  read count  the median count, in the child's reads, of the placement's tracked k-mers

Writes
  runs.tsv                      the six runs: what was audited, with which reads, the depth of one copy
  read_sets.tsv                 the nine read databases: bases, the sighting from which a k-mer is
                                kept, the peak of the count histogram
  kmers_by_parent.tsv           solid read k-mers the release lacks, by the parents' reads, before
                                any read is placed (the run's own tally)
  kmers_by_count_and_parent.tsv the same k-mers and the less frequent ones, by their count in the
                                child's reads (read from each run's scan with the package under src/)
  clusters_by_parent.tsv        placements by kind, strength and the parents' reads: how many, their
                                read counts, their shape, how close they lie to one another
  alleles_by_read_count.tsv     strong missing-allele placements by read count; beside them the
                                placements that planted errors made
  read_counts_against_planted.tsv   how many placements in neither parent's reads lie below the
                                read counts of the placements that planted errors made
  alleles_held.tsv              how many strong placements the truth genome or other assemblies of
                                the same person hold
  alleles_assemblies_agree.tsv  by how many of those assemblies each placement is held
  issues_by_parent.tsv          the same split in the issue lists
  planted_control.tsv           do errors planted in a release make calls in neither parent's reads?
  truth_genome_alleles.tsv      the strong placements on the truth genome whose k-mers are in neither
                                parent's reads
  checks.tsv                    checks of the definitions above
"""

import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from copyproof.missing import low_complexity  # noqa: E402
from copyproof.track import Scan  # noqa: E402

DATA, out = Path(sys.argv[1]), Path(sys.argv[2])
W = DATA / "work"
out.mkdir(parents=True, exist_ok=True)

# key, what is audited, run directory under DATA/work, other assemblies of the person that rode
# along as label sets, label of the truth genome where it rode along
MAIN = [
    ("HG002_hprc_year2", "HG002: the assembly hprc_year2", "rung6/hg002_calibration/hprc_year2",
     ["v2_trio", "hprc_year1", "verkko_asm_thic", "hprc_year2_v2_thic"], "truth"),
    ("HG002_truth_genome", "HG002: the truth genome itself", "hg002_v1.1/run35_trio", [], None),
    ("HG01943_release", "HG01943: HPRC release 2", "hprc_r2/HG01943", [], None),
    ("NA19240_release", "NA19240: HPRC release 2", "rung6/NA19240_short", ["verkko", "hifiasm"], None),
]
# the same genome with planted errors, run with the same reads: run directory and list of edits
PLANTED = {
    "HG002_truth_genome": ("plant/run35_trio", "plant/plant1.edits.tsv"),
    "NA19240_release": ("rung6/plant/run", "rung6/plant/plant1.edits.tsv"),
}
PAIR = {"HG002_hprc_year2": ("v2_trio", "hprc_year1"), "NA19240_release": ("verkko", "hifiasm")}  # the judges of the held-out record
SAME_READS_AS = {"HG002_hprc_year2": "HG002_truth_genome"}  # whose planted run was made with the same reads
COUNT_BINS = [(0, 8, "under 8"), (8, 9, "8 or 8.5"), (9, 11, "9 to 10.5"), (11, 13, "11 to 12.5"), (13, 17, "13 to 16.5"),
              (17, 25, "17 to 24.5"), (25, math.inf, "25 or more")]
COUNT_GROUPS = [(0, math.inf, "all"), (0, 9, "under 9"), (9, 13, "9 to 12.5"), (13, math.inf, "13 or more")]
KMER_BINS = [(3, 4, "3"), (4, 5, "4"), (5, 6, "5"), (6, 7, "6"), (7, 8, "7"), (8, 9, "8"), (9, 10, "9"), (10, 11, "10"),
             (11, 13, "11 or 12"), (13, 17, "13 to 16"), (17, 25, "17 to 24"), (25, math.inf, "25 or more")]
SIGNATURE = ["contig", "kind", "n_reads_locus", "tracked_kmers", "median_count", "reads_fwd", "reads_rev", "frac_in_mother",
             "frac_in_father", "size"]
ANCHORS = ("kmerdb", "work", "hprc_r2", "hg002", "hgsvc3", "nygc")  # first-level folders of the data root


def under_data(path: str) -> Path:
    """A path written by a run on the machine it ran on, found again under DATA."""
    parts = Path(path).parts
    at = [i for i, p in enumerate(parts) if p in ANCHORS]
    return DATA.joinpath(*parts[min(at):]) if at else Path(path)


def shown(path) -> str:
    return "DATA/" + str(under_data(str(path)).relative_to(DATA))


def load(rel: str):
    run = W / rel
    c = pd.read_csv(run / "missing_clusters.tsv", sep="\t")
    s = json.loads((run / "summary.json").read_text())
    j = json.loads((run / "run.json").read_text())
    assert c[["frac_in_mother", "frac_in_father"]].notna().all().all(), rel
    c["strength"] = np.where((c["n_reads_locus"] >= 6) & (c["tracked_kmers"] >= 5), "strong", "weak")
    c["parent_max"] = c[["frac_in_mother", "frac_in_father"]].max(axis=1)
    c["in_a_parents_reads"] = np.where(c["parent_max"] >= 0.5, "yes", "no")
    return run, c, s, j


def strong_alleles(c: pd.DataFrame) -> pd.DataFrame:
    return c[(c["kind"] == "MISSING_ALLELE") & (c["strength"] == "strong")]


def near_another(x: pd.DataFrame, within: int = 10_000) -> int:
    """How many placements have another placement of the same table within `within` bases."""
    if len(x) < 2:
        return 0
    x = x.sort_values(["contig", "start"])
    pos, con = x["start"].values.astype(float), x["contig"].values
    gap = np.where(con[1:] == con[:-1], pos[1:] - pos[:-1], np.inf)
    d = np.full(len(x), np.inf)
    d[1:] = np.minimum(d[1:], gap)
    d[:-1] = np.minimum(d[:-1], gap)
    return int((d <= within).sum())


def not_in(base: pd.DataFrame, other: pd.DataFrame):
    """Rows of `other` without a counterpart in `base` (same sequence, kind and numbers), and how
    many rows of `base` have no counterpart in `other`. Positions are not compared: an edit moves
    everything behind it."""
    left = Counter(base[SIGNATURE].itertuples(index=False, name=None))
    new = []
    for row in other[SIGNATURE].itertuples(index=False, name=None):
        if left[row] > 0:
            left[row] -= 1
            new.append(False)
        else:
            new.append(True)
    return other[np.array(new, dtype=bool)], int(sum(left.values()))


def share(k: int, n: int):
    return round(k / n, 3) if n else ""


def quartiles(v: pd.Series, name: str) -> dict:
    q = v.quantile([0.25, 0.5, 0.75])
    return {f"{name}_quartile_1": float(q[0.25]), f"{name}_median": float(q[0.5]), f"{name}_quartile_3": float(q[0.75])}


loaded = {key: load(rel) for key, _, rel, _, _ in MAIN}
planted = {key: load(rel) for key, (rel, _) in PLANTED.items()}

# ---- the runs and their reads ------------------------------------------------------------------
runs, read_sets, kmers, checks = [], [], [], []
seen_db = set()
todo = [(key, what, rel, "audited", loaded[key]) for key, what, rel, _, _ in MAIN]
todo += [(key + "_planted", dict((m[0], m[1]) for m in MAIN)[key] + ", with planted errors", rel, "control", planted[key])
         for key, (rel, _) in PLANTED.items()]
for key, what, rel, role, (run, c, s, j) in todo:
    a = c[c["kind"] == "MISSING_ALLELE"]
    db = j["databases"]
    runs.append({
        "run": key, "what": what, "role": role, "directory": f"DATA/work/{rel}", "assembly": ";".join(Path(x["path"]).name for x in j["asm"]),
        "reads": shown(db["sample"]), "mother": shown(db["mother"]), "father": shown(db["father"]),
        "other_label_sets": ";".join(k for k in db if k not in ("sample", "mother", "father")), "copyproof": j["copyproof"],
        "finished": j["finished"], "depth_per_copy": s["depth_per_copy"], "threshold": s["missing"]["threshold"],
        "placements": len(c), "missing_allele_placements": len(a), "strong_missing_allele_placements": len(strong_alleles(c)),
        "issues": len(pd.read_csv(run / "issues.tsv", sep="\t")),
    })
    sa = strong_alleles(c)
    by_event = a.groupby("twin_group")
    checks.append({
        "run": key,
        "mother_is_parent1_and_father_parent2": bool((c["frac_in_mother"] == c["frac_in_parent1"]).all()
                                                     and (c["frac_in_father"] == c["frac_in_parent2"]).all()),
        "lowest_read_count_of_a_placement": float(c["median_count"].min()),
        "allele_events_with_placements_on_both_sides_of_the_parent_cut": int((by_event["in_a_parents_reads"].nunique() > 1).sum()),
        "allele_events_with_a_strong_and_a_weak_placement": int((by_event["strength"].nunique() > 1).sum()),
        "strong_allele_placements_with_some_but_not_all_kmers_in_a_parent": int(((sa["parent_max"] > 0) & (sa["parent_max"] < 1)).sum()),
        "strong_allele_placements_with_a_parent_fraction_from_0.25_to_0.75": int(((sa["parent_max"] >= 0.25) & (sa["parent_max"] <= 0.75)).sum()),
    })
    if role != "audited":
        continue
    m = s["missing"]
    kmers.append({
        "run": key, "threshold": m["threshold"], "kmers": m["candidates"], "in_both_parents": m["candidates_in_both_parents"],
        "in_the_mother_only": m["candidates_in_parent1_only"], "in_the_father_only": m["candidates_in_parent2_only"],
        "in_neither_parent": m["candidates_in_neither_parent"],
        "share_in_neither": round(m["candidates_in_neither_parent"] / m["candidates"], 4),
    })
    for who in ("sample", "mother", "father"):
        f = under_data(db[who])
        if f in seen_db:
            continue
        seen_db.add(f)
        head = json.loads(f.with_suffix(".json").read_text())["header"]
        h = pd.read_csv(f.with_suffix(".hist.tsv"), sep="\t")
        h = h[h["count"] >= 6]  # below that the histogram is still falling from the k-mers of read errors
        peak = int(h.loc[h["n_kmers"].idxmax(), "count"])
        read_sets.append({
            "family": key.split("_")[0], "who": "child" if who == "sample" else who, "database": shown(f), "reads": head["records"],
            "bases_gb": round(head["bases"] / 1e9, 1), "kept_from_sighting": head["min_count"], "count_at_the_peak": peak,
            "half_of_it": peak / 2, "depth_per_copy_of_the_run": s["depth_per_copy"] if who == "sample" else "",
            # a k-mer present once in the genome, at a depth of half the peak, seen fewer than twice: Poisson
            "chance_of_fewer_than_2_sightings": f"{math.exp(-peak / 2) * (1 + peak / 2):.1e}",
        })
pd.DataFrame(runs).to_csv(out / "runs.tsv", sep="\t", index=False)
pd.DataFrame(read_sets).to_csv(out / "read_sets.tsv", sep="\t", index=False)
pd.DataFrame(kmers).to_csv(out / "kmers_by_parent.tsv", sep="\t", index=False)

# ---- read k-mers the release lacks, by their count in the child's reads ---------------------------
# Tracked k-mers (each stands for `sample` k-mers), without those of short tandem repeats, as the
# channel leaves them out. The channel places k-mers from the threshold up to 3.5 times the depth of
# one copy; the less frequent ones are shown beside them.
rows = []
for key, _, _, _, _ in MAIN:
    run, _, s, j = loaded[key]
    sc = Scan(str(run / "scan"))
    assert sc.reads[:3] == ["sample", "mother", "father"], sc.reads
    km, cnt = sc.novel()
    keep = ~low_complexity(km, sc.k)
    child, mo, fa = cnt[:, 0].astype(np.int64), cnt[:, 1] > 0, cnt[:, 2] > 0
    thr = s["missing"]["threshold"]
    for lo, hi, label in KMER_BINS:
        sel = keep & (child >= lo) & (child < hi)
        n0, n1, n2 = int((sel & ~mo & ~fa).sum()), int((sel & (mo ^ fa)).sum()), int((sel & mo & fa).sum())
        rows.append({"run": key, "count_in_the_childs_reads": label, "placed_by_the_channel": "yes" if lo >= thr else "no",
                     "in_neither_parent": n0, "in_one_parent": n1, "in_both_parents": n2, "share_in_neither": share(n0, n0 + n1 + n2)})
    lam = json.loads((run / "model.json").read_text())["lam"]
    cand = keep & (child >= thr) & (child <= j["parameters"]["missing"]["max_count_frac"] * lam)
    next(r for r in checks if r["run"] == key)["kmers_in_neither_parent_as_the_run_tallied"] = bool(
        int((cand & ~mo & ~fa).sum()) * sc.sample == s["missing"]["candidates_in_neither_parent"])
pd.DataFrame(rows).to_csv(out / "kmers_by_count_and_parent.tsv", sep="\t", index=False)

# ---- placements by kind, strength and the parents' reads -------------------------------------------
rows = []
for key, _, _, _, _ in MAIN:
    _, c, s, _ = loaded[key]
    lam = s["depth_per_copy"]
    for (kind, strength), g in c.groupby(["kind", "strength"]):
        for par in ("yes", "no"):
            x = g[g["in_a_parents_reads"] == par]
            if not len(x):
                continue
            row = {"run": key, "kind": kind, "strength": strength, "in_a_parents_reads": par, "placements": len(x),
                   "share_of_kind_and_strength": share(len(x), len(g)), "events": int(x["twin_group"].nunique())}
            row.update(quartiles(x["median_count"], "read_count"))
            row.update({
                "median_read_count_over_depth_per_copy": round(float(x["median_count"].median()) / lam, 2),
                "no_kmer_in_either_parent": int((x["parent_max"] == 0).sum()), "every_kmer_in_a_parent": int((x["parent_max"] == 1).sum()),
                "same_length_as_the_release": int((x["size"] == 0).sum()) if kind == "MISSING_ALLELE" else "",
            })
            row.update(quartiles(x["tracked_kmers"], "tracked_kmers"))
            row.update({
                "reads_median": float(x["n_reads_locus"].median()),
                "share_of_reads_on_the_rarer_strand_median": round(float((x[["reads_fwd", "reads_rev"]].min(axis=1) / x["n_reads_locus"]).median()), 2),
                "sequences": int(x["contig"].nunique()), "within_10_kb_of_another": near_another(x),
                "most_in_one_megabase": int(x.groupby(["contig", x["start"] // 1_000_000]).size().max()),
            })
            rows.append(row)
by_parent = pd.DataFrame(rows)
by_parent.to_csv(out / "clusters_by_parent.tsv", sep="\t", index=False)

# ---- do planted errors make calls in neither parent's reads? -----------------------------------------
control, made = [], {}
for key, (rel, edits) in PLANTED.items():
    base, pl = strong_alleles(loaded[key][1]), strong_alleles(planted[key][1])
    new, lost = not_in(base, pl)
    made[key] = new
    e = pd.read_csv(W / edits, sep="\t")
    n, k = len(new), int((new["in_a_parents_reads"] == "no").sum())
    row = {"run": key, "planted_run": key + "_planted", "edits": len(e), "edits_by_kind": ";".join(f"{a}={b}" for a, b in e["kind"].value_counts().sort_index().items())}
    for par, name in (("yes", "in_a_parent"), ("no", "in_neither")):
        row[f"unedited_{name}"] = int((base["in_a_parents_reads"] == par).sum())
        row[f"planted_{name}"] = int((pl["in_a_parents_reads"] == par).sum())
        row[f"made_by_the_edits_{name}"] = int((new["in_a_parents_reads"] == par).sum())
    row.update({"unedited_placements_not_found_again": lost, "made_by_the_edits": n, "events_made_by_the_edits": int(new["twin_group"].nunique()),
                # one-sided 95% upper bound (Clopper-Pearson) on the share in neither parent's reads
                "share_in_neither_upper_95": round(float(beta.ppf(0.95, k + 1, n - k)), 5)})
    control.append(row)
control = pd.DataFrame(control)
both = {"run": "both", "planted_run": "", "edits": int(control["edits"].sum()), "edits_by_kind": ""}
for col in control.columns[4:-1]:
    both[col] = int(control[col].sum())
both["share_in_neither_upper_95"] = round(float(beta.ppf(0.95, both["made_by_the_edits_in_neither"] + 1,
                                                         both["made_by_the_edits"] - both["made_by_the_edits_in_neither"])), 5)
pd.concat([control, pd.DataFrame([both])]).to_csv(out / "planted_control.tsv", sep="\t", index=False)

# ---- strong missing-allele placements by read count -----------------------------------------------------
rows, against = [], []
for key, _, _, _, _ in MAIN:
    s = strong_alleles(loaded[key][1])
    groups = [("in a parent's reads", s[s["in_a_parents_reads"] == "yes"]), ("in neither parent's reads", s[s["in_a_parents_reads"] == "no"])]
    if key in made:
        groups.append(("made by planted errors", made[key]))
    for name, x in groups:
        done = 0
        for lo, hi, label in COUNT_BINS:
            n = int(((x["median_count"] >= lo) & (x["median_count"] < hi)).sum())
            done += n
            rows.append({"run": key, "placements_are": name, "of": len(x), "read_count": label, "placements": n, "share": share(n, len(x)),
                         "share_up_to_here": share(done, len(x))})
    # against the read counts of placements that planted errors made (alleles that the genome holds
    # and the edited release lacks), where the same reads were used for a planted run
    ref_key = SAME_READS_AS.get(key, key)
    neither = groups[1][1]
    refs = [("this run's placements in a parent's reads", groups[0][1])]
    if ref_key in made:
        refs.append((f"placements made by planted errors ({ref_key}_planted)", made[ref_key]))
    for ref_name, ref in refs:
        for q, q_name in ((0.05, "its 5th percentile"), (0.25, "its lower quartile"), (0.5, "its median")):
            cut = float(ref["median_count"].quantile(q))
            below = int((neither["median_count"] < cut).sum())
            against.append({"run": key, "reference": ref_name, "reference_placements": len(ref), "cut": q_name, "read_count_at_the_cut": cut,
                            "reference_below_the_cut": int((ref["median_count"] < cut).sum()),
                            "in_neither_parents_reads": len(neither), "of_them_below_the_cut": below, "share": share(below, len(neither))})
pd.DataFrame(rows).to_csv(out / "alleles_by_read_count.tsv", sep="\t", index=False)
pd.DataFrame(against).to_csv(out / "read_counts_against_planted.tsv", sep="\t", index=False)

# ---- what other assemblies of the same person hold --------------------------------------------------------
rows, agree = [], []
for key, _, _, judges, truth in MAIN:
    labels = ([truth] if truth else []) + judges
    if not labels:
        continue
    s = strong_alleles(loaded[key][1])
    held = {lab: (s[f"frac_in_{lab}"] >= 0.5) for lab in labels}
    combos = dict(held)
    a, b = PAIR[key]
    combos[f"{a} or {b}"] = held[a] | held[b]
    combos[f"{a} and {b}"] = held[a] & held[b]
    part = int(sum(((s[f"frac_in_{lab}"] > 0) & (s[f"frac_in_{lab}"] < 1)).sum() for lab in labels))
    checks_row = next(r for r in checks if r["run"] == key)
    checks_row["strong_allele_placements_held_in_part_by_an_assembly"] = part
    n_hold = sum(h.astype(int) for h in held.values())
    for par in ("yes", "no", "yes or no"):
        of_class = (s["in_a_parents_reads"] == par) | (par == "yes or no")
        for lo, hi, label in COUNT_GROUPS:
            sel = of_class & (s["median_count"] >= lo) & (s["median_count"] < hi)
            for lab, h in combos.items():
                rows.append({"run": key, "in_a_parents_reads": par, "read_count": label, "held_by": lab, "placements": int(sel.sum()),
                             "held": int((sel & h).sum()), "share": share(int((sel & h).sum()), int(sel.sum()))})
        if par == "yes or no":
            continue
        sel = s["in_a_parents_reads"] == par
        for k in range(len(labels) + 1):
            agree.append({"run": key, "in_a_parents_reads": par, "assemblies": ";".join(labels), "held_by_this_many": k,
                          "placements": int((sel & (n_hold == k)).sum()), "share": share(int((sel & (n_hold == k)).sum()), int(sel.sum()))})
pd.DataFrame(rows).to_csv(out / "alleles_held.tsv", sep="\t", index=False)
pd.DataFrame(agree).to_csv(out / "alleles_assemblies_agree.tsv", sep="\t", index=False)
pd.DataFrame(checks).to_csv(out / "checks.tsv", sep="\t", index=False)

# ---- the issue lists ------------------------------------------------------------------------------------------
rows = []
for key, _, _, _, _ in MAIN:
    iss = pd.read_csv(loaded[key][0] / "issues.tsv", sep="\t")
    m = iss[iss["class"].str.startswith("MISSING_")].copy()
    p1 = m["evidence"].str.extract(r"in_parent1=([0-9.]+)")[0].astype(float)
    p2 = m["evidence"].str.extract(r"in_parent2=([0-9.]+)")[0].astype(float)
    assert p1.notna().all() and p2.notna().all(), key
    # for MISSING_HAPLOTYPE the evidence gives the mean over the alleles of the run
    m["in_a_parents_reads"] = np.where((p1 >= 0.5) | (p2 >= 0.5), "yes", "no")
    g = m.groupby(["class", "role", "confidence", "in_a_parents_reads"]).size().unstack("in_a_parents_reads", fill_value=0)
    for (cls, role, conf), r in g.iterrows():
        yes, no = int(r.get("yes", 0)), int(r.get("no", 0))
        rows.append({"run": key, "class": cls, "role": role, "confidence": conf, "issues": yes + no, "in_a_parents_reads": yes,
                     "in_neither_parents_reads": no, "share_in_neither": share(no, yes + no)})
    for cls, conf, sel in (("MISSING_ALLELE", "high", (m["class"] == "MISSING_ALLELE") & (m["confidence"] == "high")),
                           ("every MISSING_ class", "high or low", pd.Series(True, index=m.index)),
                           ("the whole issue list", "high or low", None)):
        no = int((m["in_a_parents_reads"] == "no").sum()) if sel is None else int((sel & (m["in_a_parents_reads"] == "no")).sum())
        n = len(iss) if sel is None else int(sel.sum())
        rows.append({"run": key, "class": cls, "role": "locus or twin", "confidence": conf, "issues": n,
                     "in_a_parents_reads": "" if sel is None else n - no, "in_neither_parents_reads": no, "share_in_neither": share(no, n)})
issues = pd.DataFrame(rows)
issues.to_csv(out / "issues_by_parent.tsv", sep="\t", index=False)

# ---- the placements on the truth genome ---------------------------------------------------------------------------
t = strong_alleles(loaded["HG002_truth_genome"][1])
t = t[t["in_a_parents_reads"] == "no"].sort_values(["twin_group", "contig"])
t[["twin_group", "contig", "start", "end", "size", "n_reads_locus", "reads_fwd", "reads_rev", "tracked_kmers", "median_count",
   "frac_in_mother", "frac_in_father", "other_loci"]].to_csv(out / "truth_genome_alleles.tsv", sep="\t", index=False)

# ---- to the terminal ------------------------------------------------------------------------------------------------
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
sa = by_parent[(by_parent["kind"] == "MISSING_ALLELE") & (by_parent["strength"] == "strong")]
print(sa[["run", "in_a_parents_reads", "placements", "share_of_kind_and_strength", "events", "read_count_median",
          "no_kmer_in_either_parent", "same_length_as_the_release", "within_10_kb_of_another"]].to_string(index=False))
print(pd.read_csv(out / "planted_control.tsv", sep="\t")[["run", "made_by_the_edits", "made_by_the_edits_in_a_parent", "made_by_the_edits_in_neither",
                                                              "unedited_placements_not_found_again", "share_in_neither_upper_95"]].to_string(index=False))
print(issues[issues["role"] == "locus or twin"].to_string(index=False))
