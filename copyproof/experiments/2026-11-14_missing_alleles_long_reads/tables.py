#!/usr/bin/env python3
"""The tables of this record: what the strong missing-allele placements are, read from the long
reads each release was built from and from the changed bases themselves.

  tables.py DATA OUT_DIR [RUN ...]

With run keys after OUT_DIR only those runs are tabulated (for a trial; the record's tables are made
without).

DATA is the data root. Read: four finished runs under DATA/work (cluster table, `missing/` folder,
scan side files, `run.json`, `model.json`), the releases they audited (through `cpk extract`), and
the long reads that `fetch_long_reads.py` put under DATA/work/missing_alleles_long_reads. Written
under DATA, all in DATA/work/missing_alleles_long_reads: <run>/release.regions.* (the release's
sequence around each placement, kept so that a second run does not load the assemblies again),
`placements.tsv` (one row per strong missing-allele placement, with its change and its long reads;
1.4 MB, so not in the record) and `other_rows_of_the_events.tsv` (the long reads at the other rows of
the same events). No network is used.

Terms as in the record of 2026-11-14 (`../2026-11-14_missing_alleles_by_parent`): placement, event,
strong, in a parent's reads, read count. How a placement is sorted from its long reads was fixed
in NOTEBOOK entry u before any long read was searched; `placements.py` holds that rule as code.

Writes
  runs.tsv                        the four runs and, for the two releases, what was fetched
  checks.tsv                      every row of each cluster table against what was recovered for it;
                                  the changed bases against a second way of deriving them
  long_reads_by_parent.tsv        placements and events per long-read group, by the parents' reads
  long_reads_events.tsv           events sorted from all their placements
  long_reads_event_patterns.tsv   which groups the placements of one event fall into together
  long_reads_by_read_count.tsv    the groups by the read count of the short reads
  long_reads_by_other_assemblies.tsv   NA19240: the groups by what the HGSVC3 assemblies hold
  long_reads_events_by.tsv        the same two splits with each event counted once
  long_reads_share_steps.tsv      the share of long reads with the allele in finer steps
  long_reads_reads.tsv            the long reads themselves, summed over placements
  long_reads_aligned_base.tsv     one changed base: the k-mers of a read against its aligned base
  long_reads_by_aligned_base.tsv  one changed base: the groups by k-mers against the groups by aligned base
  share_of_all_reads.tsv          events: the allele's share of all long reads at both places, beside the
                                  read count of the short reads over two copies
  long_reads_by_change.tsv        events: the groups by the kind of change
  long_reads_at_the_nearer_placement.tsv   events with one placement of one changed base and one of a larger
                                  difference: at which the long reads hold the allele
  changed_bases.tsv               how many placements are one changed base, and what the others are
  changed_bases_detail.tsv        the same with every reason as the script words it
  size_column.tsv                 what the cluster table's `size` of 0 covers
  substitutions.tsv               the six kinds of substitution
  substitutions_with_neighbours.tsv   the six kinds with the base before and the base after
  substitutions_by_long_reads.tsv the six kinds by long-read group
  substitutions_beside_the_cluster_table.tsv   the six kinds by read count, strands and what another assembly holds
"""

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import placements as P  # noqa: E402

DATA, out = Path(sys.argv[1]), Path(sys.argv[2])
W = DATA / "work"
out.mkdir(parents=True, exist_ok=True)
CLASS = {"yes": "in a parent's reads", "no": "in neither parent's reads"}
KINDS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CHANGE_GROUPS = ["one changed base", "two or more changed bases less than 31 apart", "another length: 1 base", "another length: 2 to 30 bases",
                 "another length: more than 30 bases", "not determined"]


def share(k: int, n: int):
    return round(k / n, 3) if n else ""


def shown(path) -> str:
    return "DATA/" + str(P.under_data(DATA, str(path)).relative_to(DATA))


def change_group(resolved: str) -> str:
    """The reasons of placements.changed_base, in six groups."""
    if resolved == "one changed base":
        return resolved
    if "changed bases within" in resolved and "more than" not in resolved:
        return "two or more changed bases less than 31 apart"
    if "more than the release" in resolved or "fewer than the release" in resolved:
        n = int(resolved.split("have ")[1].split(" ")[0])
        return "another length: 1 base" if n == 1 else ("another length: 2 to 30 bases" if n <= 30 else "another length: more than 30 bases")
    return "not determined"


def held_by(verkko: float, hifiasm: float) -> str:
    v, h = verkko >= 0.5, hifiasm >= 0.5
    return "both" if v and h else ("the Verkko assembly only" if v else ("the hifiasm assembly only" if h else "neither"))


# ---- every placement: recover, derive the change, count the long reads --------------------------------
runs, checks, rows, mates_rows = [], [], [], []
ONLY = sys.argv[3:]
for key, what, rel, script in P.RUNS:
    if ONLY and key not in ONLY:
        continue
    run = W / rel
    table = P.load_table(run)
    rec = P.recover(run)
    chk = P.compare(table, rec)
    fields = [c for c in chk.columns if c not in ("cluster", "kind", "strength")]
    sa = P.strong_alleles(table)
    events = set(sa["twin_group"])
    mates = table[table["twin_group"].isin(events) & ~table["cluster"].isin(sa["cluster"])]
    by_event = sa.groupby("twin_group")["in_a_parents_reads"].nunique()
    todo = pd.concat([sa, mates])
    lengths = dict(pd.read_csv(run / "scan.contigs.tsv", sep="\t", dtype={"name": str})[["name", "length"]].itertuples(index=False))
    regions = pd.DataFrame({"name": todo["cluster"].astype(str), "contig": todo["contig"], "start": (todo["start"] - P.FLANK).clip(lower=0),
                            "end": [min(int(e) + P.FLANK, int(lengths[c])) for e, c in zip(todo["end"], todo["contig"])]})
    seqs = P.release_sequences(DATA, run, regions, W / P.WORK / key / "release")
    wstart = dict(zip(regions["name"], regions["start"]))
    j = json.loads((run / "run.json").read_text())
    row = {"run": key, "what": what, "directory": f"DATA/work/{rel}", "assembly": ";".join(Path(a["path"]).name for a in j["asm"]),
           "reads": shown(j["databases"]["sample"]), "copyproof": j["copyproof"], "finished": j["finished"],
           "strong_missing_allele_placements": len(sa), "in_a_parents_reads": int((sa["in_a_parents_reads"] == "yes").sum()),
           "in_neither_parents_reads": int((sa["in_a_parents_reads"] == "no").sum()), "events": len(events),
           "events_in_a_parents_reads": int(sa.loc[sa["in_a_parents_reads"] == "yes", "twin_group"].nunique()),
           "events_in_neither_parents_reads": int(sa.loc[sa["in_a_parents_reads"] == "no", "twin_group"].nunique()),
           "other_rows_of_those_events": len(mates)}
    bam = None
    if script:
        folder = W / P.WORK / key
        bam = folder / "long_reads.bam"
        f = json.loads((folder / "fetch.json").read_text())
        # how the public file says it was made: the aligner lines of its header
        header = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True).stdout.splitlines()
        aligned = [ln for ln in header if ln.startswith("@PG") and "\tPN:minimap2" in ln]
        cl = aligned[0].split("\tCL:")[1].split(" ")
        row.update({"long_reads": f["address"], "aligner": " ".join(x for x in cl[:-2] if not x.startswith("-t")) + f" (version {aligned[0].split(chr(9) + 'VN:')[1].split(chr(9))[0]})",
                    "aligned_to": cl[-2], "read_files_aligned": len(aligned), "first_read_file": Path(cl[-1]).name,
                    # the alignment was made to both haplotypes at once when its header lists every sequence of the release
                    "sequences_of_the_release": len(lengths),
                    "of_them_in_the_header_of_the_alignment": sum(1 for ln in header if ln.startswith("@SQ") and ln.split("\tSN:")[1].split("\t")[0] in lengths),
                    "sequences_in_the_header": sum(1 for ln in header if ln.startswith("@SQ")),
                    "fetched": f["started"], "fetch_seconds": f["seconds"], "regions_asked_for": f["regions"],
                    "records_fetched": f["records"], "of_them_secondary_or_supplementary": f["records"] - f["records_not_secondary_or_supplementary"],
                    "read_bases_fetched_mb": round(f["read_bases"] / 1e6, 1), "bam_mb": round(f["bam_bytes"] / 1e6, 1), "kept_in": f"DATA/work/{P.WORK}/{key}"})
    row["depth_per_copy"] = round(rec.lam, 2)
    runs.append(row)

    second = Counter()
    for r in todo.itertuples(index=False):
        c = int(r.cluster)
        window, w0 = seqs[str(c)], int(wstart[str(c)])
        strong = r.kind == "MISSING_ALLELE" and r.strength == "strong"
        one = {"run": key, "cluster": c, "sequence": r.contig, "start": int(r.start), "end": int(r.end), "size": int(r.size),
               "event": int(r.twin_group), "in_a_parents_reads": r.in_a_parents_reads, "reads": int(r.n_reads_locus),
               "tracked_kmers": int(r.tracked_kmers), "read_count": float(r.median_count)}
        if not strong:
            one.update({"kind": r.kind, "strength": r.strength})
        cb = None
        if strong:
            one["read_count_group"] = P.read_count_group(float(r.median_count))
            one["reads_on_the_rarer_strand"] = int(min(r.reads_fwd, r.reads_rev))
            if "frac_in_verkko" in table.columns:
                one["held_by_the_hgsvc3_assemblies"] = held_by(float(r.frac_in_verkko), float(r.frac_in_hifiasm))
                one["held_elsewhere"] = "yes" if max(float(r.frac_in_verkko), float(r.frac_in_hifiasm)) >= 0.5 else "no"
            elif "frac_in_truth" in table.columns and key != "HG002_truth_genome":
                one["held_elsewhere"] = "yes" if float(r.frac_in_truth) >= 0.5 else "no"
            cb = P.changed_base(rec, c, int(r.start), int(r.end), window, w0)
            one["change"] = cb["resolved"]
            one["change_group"] = change_group(cb["resolved"])
            for col in ("release_position", "release_has", "reads_have", "before", "after"):
                one[col] = cb[col]
            if cb["resolved"] == "one changed base":
                one["substitution"], one["base_before"], one["base_after"] = P.substitution(cb["release_base"], cb["reads_base"], cb["before"], cb["after"])
            # the second way: the reads laid over one another, compared with the release between two anchors
            ref_mid, alt_mid, at = P.between_anchors(rec, *rec.members[c][:1], rec.members[c][2], int(r.start), window, w0)
            if cb["resolved"] == "one changed base":
                second["one changed base by the k-mers"] += 1
                i = cb["release_position"] - w0 - at if at is not None else -1
                if at is None:
                    second["no sequence"] += 1
                elif (ref_mid, alt_mid, i) == (cb["release_base"], cb["reads_base"], 0):
                    second["the same base, alone"] += 1
                elif len(ref_mid) == len(alt_mid) and 0 <= i < len(alt_mid) and (ref_mid[i], alt_mid[i]) == (cb["release_base"], cb["reads_base"]):
                    second["the same base, beside others"] += 1  # the anchors lie further out than the placement's k-mers reach
                elif len(ref_mid) != len(alt_mid):
                    second["another length"] += 1
                else:
                    second["another base"] += 1
        if bam is not None:
            allele = P.decode_kmers(rec.kmers[rec.members[c][1]], rec.k)
            release = P.release_kmers_there(rec, c, int(r.start), int(r.end), window, w0)
            at = cb["release_position"] if cb is not None and cb["resolved"] == "one changed base" else -1
            lr = P.long_reads_at(bam, r.contig, int(r.start), int(r.end), allele, release, rec.k, changed_at=at)
            either, sh = P.share_of(lr["hold_the_allele"], lr["hold_the_release"], lr["hold_both"])
            one["release_kmers_there"] = len(release)
            one.update({f"long_reads_{k}": v for k, v in lr.items() if k != "bases"})
            # where the release has no k-mer of that kind, no read can hold its sequence: too few to say
            one.update({"long_reads_hold_either": either, "share_with_the_allele": round(sh, 3) if either and release else "",
                        "long_read_group": P.group_of(lr["hold_the_allele"], either) if release else P.GROUPS[4], "aligned_bases": lr["bases"]})
            if cb is not None and cb["resolved"] == "one changed base":
                b = Counter()
                for item in filter(None, lr["bases"].split(";")):
                    base, rest = item.split(":")
                    b["reads" if base == cb["reads_base"] else ("release" if base == cb["release_base"] else "other")] += int(rest.split("=")[1])
                one.update({"aligned_base_is_the_reads": b["reads"], "aligned_base_is_the_releases": b["release"], "aligned_base_is_another": b["other"],
                            "group_by_aligned_base": P.group_of(b["reads"], b["reads"] + b["release"])})
        (rows if strong else mates_rows).append(one)

    checks.append({
        "run": key, "rows_of_the_cluster_table": len(table), "rows_recovered_under_the_same_cluster_number": int(chk["cluster"].isin(list(rec.members)).sum()),
        **{f"same_{f}": int(chk[f].sum()) for f in fields},
        "strong_missing_allele_rows": len(sa), **{f"strong_same_{f}": int(chk.loc[chk["cluster"].isin(sa["cluster"]), f].sum()) for f in ("reads", "tracked_kmers", "read_count")},
        "kmer_list_of_the_missing_folder_is_the_scans": rec.kmer_list_is_the_runs,
        "strong_allele_events_on_both_sides_of_the_parent_cut": int((by_event > 1).sum()),
        "one_changed_base_by_the_kmers": second["one changed base by the k-mers"],
        "second_way_the_same_base_alone": second["the same base, alone"],
        "second_way_the_same_base_beside_other_differences_between_the_anchors": second["the same base, beside others"],
        "second_way_another_length_between_the_anchors": second["another length"],
        "second_way_another_base": second["another base"],
        "second_way_no_sequence_anchored": second["no sequence"],
    })
    print(f"{key}: {len(sa)} strong missing-allele placements, {len(mates)} other rows of their events", flush=True)

S = pd.DataFrame(rows)
M = pd.DataFrame(mates_rows)
# An allele belongs to one haplotype. At the other haplotype's place it differs from the release by
# its own change and by whatever differs between the two haplotypes there. So an event is counted
# under the smallest difference among its strong placements, in the order of CHANGE_GROUPS.
nearest = S.groupby(["run", "event"])["change_group"].agg(lambda v: min(v, key=CHANGE_GROUPS.index)).rename("event_change").reset_index()
S = S.merge(nearest, on=["run", "event"])
pd.DataFrame(runs).convert_dtypes().to_csv(out / "runs.tsv", sep="\t", index=False)
pd.DataFrame(checks).to_csv(out / "checks.tsv", sep="\t", index=False)
# one row per placement is 1.4 MB: it is written beside the fetched reads, not into the record
lead = ["run", "cluster", "sequence", "start", "end", "size", "event", "in_a_parents_reads", "reads", "reads_on_the_rarer_strand", "tracked_kmers", "read_count",
        "read_count_group", "held_by_the_hgsvc3_assemblies", "held_elsewhere", "change_group", "change", "release_position", "release_has", "reads_have", "before", "after", "substitution",
        "base_before", "base_after"]
if not ONLY:
    S[lead + [c for c in S.columns if c not in lead]].to_csv(W / P.WORK / "placements.tsv", sep="\t", index=False)
M = M[M["long_read_group"].notna()].drop(columns=["aligned_bases"])  # the two releases
if not ONLY:
    M.to_csv(W / P.WORK / "other_rows_of_the_events.tsv", sep="\t", index=False)

L = S[S["long_read_group"].notna()].copy()  # the two releases
RELEASES = [k for k, _, _, script in P.RUNS if script and (not ONLY or k in ONLY)]

# ---- an event from all its rows -----------------------------------------------------------------------------
ev = []
for key in RELEASES:
    a, b = L[L["run"] == key], M[M["run"] == key]
    groups = {e: list(g) for e, g in a.groupby("event")["long_read_group"]}
    for e, g in b.groupby("event")["long_read_group"]:
        groups[e] = groups[e] + list(g)
    cls = a.groupby("event")["in_a_parents_reads"].first()
    n_strong = a.groupby("event").size()
    for e, gs in groups.items():
        ev.append({"run": key, "event": e, "in_a_parents_reads": cls[e], "strong_placements": int(n_strong[e]), "rows": len(gs),
                   "event_group": P.event_group(gs), "pattern": " + ".join(sorted(gs, key=P.GROUPS.index))})
E = pd.DataFrame(ev)
L = L.merge(E[["run", "event", "event_group"]], on=["run", "event"])


ORDER = {v: i for i, v in enumerate([r[0] for r in P.RUNS] + list(CLASS.values()) + ["under 9", "9 to 12.5", "13 or more", "both",
                                    "the Verkko assembly only", "the hifiasm assembly only", "neither"])}


def by_group(frame: pd.DataFrame, keys: list[str], extra=None) -> pd.DataFrame:
    """Placements and events per long-read group within each combination of `keys`."""
    res = []
    parts = sorted(frame.groupby(keys, sort=False), key=lambda kv: [ORDER.get(v, 99) for v in kv[0]])
    for vals, g in parts:
        for name in P.GROUPS:
            x = g[g["long_read_group"] == name]
            row = dict(zip(keys, vals))
            row.update({"long_reads": name, "placements": len(x), "of": len(g), "share": share(len(x), len(g)), "events": int(x["event"].nunique()),
                        "of_events": int(g["event"].nunique())})
            if extra:
                row.update(extra(x))
            res.append(row)
    return pd.DataFrame(res)


def reads_of(x: pd.DataFrame) -> dict:
    either, allele = int(x["long_reads_hold_either"].sum()), int(x["long_reads_hold_the_allele"].sum())
    low = int((x["long_reads_covering_with_mapping_quality_under_20"] >= 0.5 * x["long_reads_covering"].clip(lower=1)).sum())
    return {"long_reads_that_hold_either": either, "of_them_hold_the_allele": allele, "of_them_hold_the_release_only": either - allele,
            "share_of_those_reads_with_the_allele": share(allele, either),
            "median_reads_that_hold_either": float(x["long_reads_hold_either"].median()) if len(x) else "",
            # where the two haplotypes are alike over a read's length, the aligner cannot tell them apart
            "placements_where_half_the_reads_or_more_have_mapping_quality_under_20": low}


L["class"] = L["in_a_parents_reads"].map(CLASS)
order = pd.CategoricalDtype(RELEASES, ordered=True)
L["run"] = L["run"].astype(order)
L = L.sort_values(["run", "in_a_parents_reads", "cluster"], ascending=[True, False, True])
L["run"] = L["run"].astype(str)
by_parent = by_group(L, ["run", "class"], reads_of)
by_parent.to_csv(out / "long_reads_by_parent.tsv", sep="\t", index=False)
by_group(L, ["run", "class", "read_count_group"]).to_csv(out / "long_reads_by_read_count.tsv", sep="\t", index=False)
by_group(L[L["run"] == "NA19240_release"], ["run", "class", "held_by_the_hgsvc3_assemblies"]).to_csv(out / "long_reads_by_other_assemblies.tsv", sep="\t", index=False)

E["class"] = E["in_a_parents_reads"].map(CLASS)
res = []
for (key, cls), g in E.groupby(["run", "class"], sort=False):
    for name in P.GROUPS:
        n = int((g["event_group"] == name).sum())
        res.append({"run": key, "class": cls, "long_reads": name, "events": n, "of": len(g), "share": share(n, len(g))})
events = pd.DataFrame(res).sort_values(["run", "class"], key=lambda v: v.map(ORDER), kind="stable")
events.to_csv(out / "long_reads_events.tsv", sep="\t", index=False)
pat = E.groupby(["run", "class", "rows", "pattern"]).size().rename("events").reset_index().sort_values(["run", "class", "events"], ascending=[True, True, False])
pat.to_csv(out / "long_reads_event_patterns.tsv", sep="\t", index=False)

# How much of all the DNA holds the allele: an event's long reads summed over all its rows (both
# haplotypes' places), beside the read count of its short reads over two copies. Events with 20 long
# reads or more that hold the allele or the release's sequence.
lam = {r["run"]: r["depth_per_copy"] for r in runs}
cols = ["run", "event", "long_reads_hold_the_allele", "long_reads_hold_either"]
sums = pd.concat([L[cols], M[cols]]).groupby(["run", "event"]).sum().reset_index()
counts = L.groupby(["run", "event"])["read_count"].median().rename("read_count").reset_index()
V = E.merge(sums, on=["run", "event"]).merge(counts, on=["run", "event"])
res = []
for (key, cls), g in sorted(V.groupby(["run", "class"], sort=False), key=lambda kv: [ORDER[v] for v in kv[0]]):
    x = g[g["long_reads_hold_either"] >= 20]
    lr = x["long_reads_hold_the_allele"] / x["long_reads_hold_either"]
    sr = x["read_count"] / (2 * lam[key])
    row = {"run": key, "class": cls, "events": len(g), "events_with_20_long_reads_or_more": len(x),
           "long_reads_of_those_events": int(x["long_reads_hold_either"].sum()), "of_them_hold_the_allele": int(x["long_reads_hold_the_allele"].sum()),
           "share_summed": share(int(x["long_reads_hold_the_allele"].sum()), int(x["long_reads_hold_either"].sum()))}
    for q, name in ((0.25, "quartile_1"), (0.5, "median"), (0.75, "quartile_3")):
        row[f"share_of_the_events_long_reads_{name}"] = round(float(lr.quantile(q)), 3) if len(x) else ""
    for q, name in ((0.25, "quartile_1"), (0.5, "median"), (0.75, "quartile_3")):
        row[f"short_read_count_over_two_copies_{name}"] = round(float(sr.quantile(q)), 3) if len(x) else ""
    row["depth_per_copy_of_the_short_reads"] = lam[key]
    res.append(row)
pd.DataFrame(res).to_csv(out / "share_of_all_reads.tsv", sep="\t", index=False)

# events, each once, by the read count of the short reads and by what the HGSVC3 assemblies hold
first = L.groupby(["run", "event"]).agg(read_count=("read_count", "median"), held=("held_by_the_hgsvc3_assemblies", "first")).reset_index()
first["read_count_group"] = first["read_count"].map(P.read_count_group)
X = E.merge(first, on=["run", "event"])
res = []
for split, col, values in (("read count of the short reads", "read_count_group", ["under 9", "9 to 12.5", "13 or more"]),
                           ("the HGSVC3 assemblies that hold the allele", "held", ["both", "the Verkko assembly only", "the hifiasm assembly only", "neither"])):
    for (key, cls), g in sorted(X.groupby(["run", "class"], sort=False), key=lambda kv: [ORDER[v] for v in kv[0]]):
        if g[col].isna().all():
            continue
        for v in values:
            x = g[g[col] == v]
            row = {"run": key, "class": cls, "events_by": split, "value": v, "events": len(x)}
            for name in P.GROUPS:
                row[name] = int((x["event_group"] == name).sum())
            res.append(row)
pd.DataFrame(res).to_csv(out / "long_reads_events_by.tsv", sep="\t", index=False)

# the event's long reads by the event's change
E = E.merge(nearest, on=["run", "event"])
res = []
for (key, cls), g in sorted(E.groupby(["run", "class"], sort=False), key=lambda kv: [ORDER[v] for v in kv[0]]):
    for change in CHANGE_GROUPS:
        x = g[g["event_change"] == change]
        row = {"run": key, "class": cls, "change": change, "events": len(x)}
        for name in P.GROUPS:
            row[name] = int((x["event_group"] == name).sum())
        res.append(row)
pd.DataFrame(res).to_csv(out / "long_reads_by_change.tsv", sep="\t", index=False)

# Two sources that do not know of each other: where an event has one placement of one changed base and
# another placement with a larger difference, at which of the two do the long reads hold the allele?
res = []
HOLD = set(P.GROUPS[1:4])
for (key, cls), g in sorted(L.groupby(["run", "class"], sort=False), key=lambda kv: [ORDER[v] for v in kv[0]]):
    tally = Counter()
    for _, e in g.groupby("event"):
        near, far = e[e["change_group"] == CHANGE_GROUPS[0]], e[~e["change_group"].isin([CHANGE_GROUPS[0], "not determined"])]
        if len(e) != 2 or len(near) != 1 or len(far) != 1:
            continue
        a, b = near["long_read_group"].iloc[0], far["long_read_group"].iloc[0]
        if P.GROUPS[4] in (a, b):
            tally["too few reads at one of the two"] += 1
        else:
            tally[("at the placement of one changed base only" if b not in HOLD else "at both") if a in HOLD else
                  ("at the other placement only" if b in HOLD else "at neither")] += 1
    for name in ("at the placement of one changed base only", "at the other placement only", "at both", "at neither", "too few reads at one of the two"):
        res.append({"run": key, "class": cls, "events_with_one_placement_of_one_changed_base_and_one_of_a_larger_difference": sum(tally.values()),
                    "long_reads_hold_the_allele": name, "events": tally[name]})
pd.DataFrame(res).to_csv(out / "long_reads_at_the_nearer_placement.tsv", sep="\t", index=False)

# ---- the share in finer steps ---------------------------------------------------------------------------
res = []
for (key, cls), g in L.groupby(["run", "class"], sort=False):
    g = g[g["release_kmers_there"] > 0]
    has = g[g["long_reads_hold_either"] > 0]
    step = has.apply(lambda r: P.step_of(r["long_reads_hold_the_allele"] / r["long_reads_hold_either"]), axis=1) if len(has) else pd.Series(dtype=str)
    enough = has["long_reads_hold_either"] >= P.ENOUGH
    for name in P.SHARE_STEPS:
        res.append({"run": key, "class": cls, "share_of_long_reads_with_the_allele": name, "placements_with_10_reads_or_more": int(((step == name) & enough).sum()),
                    "of": int(enough.sum()), "placements_with_1_to_9_reads": int(((step == name) & ~enough).sum())})
    res.append({"run": key, "class": cls, "share_of_long_reads_with_the_allele": "no read holds either", "placements_with_10_reads_or_more": "",
                "of": "", "placements_with_1_to_9_reads": int((g["long_reads_hold_either"] == 0).sum())})
pd.DataFrame(res).to_csv(out / "long_reads_share_steps.tsv", sep="\t", index=False)

# ---- the long reads themselves ------------------------------------------------------------------------------
res = []
for (key, cls), g in L.groupby(["run", "class"], sort=False):
    row = {"run": key, "class": cls, "placements": len(g)}
    for col in ("records", "left_out_secondary_or_supplementary", "not_covering", "covering", "hold_the_allele", "hold_the_release", "hold_both", "hold_neither",
                "allele_only_outside_the_window", "hold_every_allele_kmer", "hold_under_half_of_the_allele_kmers", "covering_with_mapping_quality_under_20",
                "allele_holders_with_mapping_quality_under_20"):
        row[col] = int(g[f"long_reads_{col}"].sum())
    row.update({"covering_per_placement_quartile_1": float(g["long_reads_covering"].quantile(0.25)), "covering_per_placement_median": float(g["long_reads_covering"].median()),
                "covering_per_placement_quartile_3": float(g["long_reads_covering"].quantile(0.75)),
                "placements_without_a_release_kmer_there": int((g["release_kmers_there"] == 0).sum()),
                "placements_with_fewer_than_10_covering": int((g["long_reads_covering"] < P.ENOUGH).sum())})
    res.append(row)
pd.DataFrame(res).to_csv(out / "long_reads_reads.tsv", sep="\t", index=False)

# ---- one changed base: k-mers against the aligned base ----------------------------------------------------------
one = L[L["change_group"] == "one changed base"]
res = []
for (key, cls), g in one.groupby(["run", "class"], sort=False):
    tally = Counter()
    for r in g.itertuples(index=False):
        rel_base, reads_base = r.release_has, r.reads_have
        for item in filter(None, r.aligned_bases.split(";")):
            base, rest = item.split(":")
            what, n = rest.split("=")
            col = "the reads' base" if base == reads_base else ("the release's base" if base == rel_base else ("no base (deleted)" if base == "-" else "another base"))
            tally[(what, col)] += int(n)
    for what, label in (("allele", "holds k-mers of the allele"), ("release", "holds k-mers of the release only"), ("neither", "holds neither")):
        row = {"run": key, "class": cls, "placements": len(g), "long_read": label}
        for col in ("the reads' base", "the release's base", "another base", "no base (deleted)"):
            row[f"aligned_base_is_{col}"] = tally[(what, col)]
        res.append(row)
pd.DataFrame(res).to_csv(out / "long_reads_aligned_base.tsv", sep="\t", index=False)
res = []
for (key, cls), g in one.groupby(["run", "class"], sort=False):
    for a in P.GROUPS:
        row = {"run": key, "class": cls, "by_kmers": a, "placements": int((g["long_read_group"] == a).sum())}
        for b in P.GROUPS:
            row[f"by_aligned_base: {b}"] = int(((g["long_read_group"] == a) & (g["group_by_aligned_base"] == b)).sum())
        res.append(row)
pd.DataFrame(res).to_csv(out / "long_reads_by_aligned_base.tsv", sep="\t", index=False)

# ---- the changed bases ---------------------------------------------------------------------------------------------
S["class"] = S["in_a_parents_reads"].map(CLASS)
S["run"] = S["run"].astype(pd.CategoricalDtype([r[0] for r in P.RUNS], ordered=True))
S = S.sort_values(["run", "in_a_parents_reads", "cluster"], ascending=[True, False, True])
S["run"] = S["run"].astype(str)
res, detail, sizes = [], [], []
for (key, cls), g in S.groupby(["run", "class"], sort=False):
    once = g.drop_duplicates("event")
    for name in CHANGE_GROUPS:
        x = g[g["change_group"] == name]
        n = int((once["event_change"] == name).sum())
        res.append({"run": key, "class": cls, "change": name, "placements": len(x), "of": len(g), "share": share(len(x), len(g)),
                    "events_by_their_smallest_difference": n, "of_events": len(once), "share_of_events": share(n, len(once))})
    reason = g["change"].str.replace(r"\d+", "N", regex=True)
    for name, n in reason.value_counts().items():
        x = g[reason == name]
        detail.append({"run": key, "class": cls, "change_group": x["change_group"].iloc[0], "change_as_worded_by_the_script": name, "placements": int(n),
                       "events": int(x["event"].nunique())})
    zero = g[g["size"] == 0]
    longer = zero[zero["change_group"].str.startswith("another length")]
    d = (longer["reads_have"].str.len() - longer["release_has"].str.len()).abs()
    sizes.append({"run": key, "class": cls, "placements": len(g), "size_is_0_in_the_cluster_table": len(zero),
                  "of_them_one_changed_base": int((zero["change_group"] == "one changed base").sum()),
                  "of_them_two_or_more_changed_bases": int((zero["change_group"] == CHANGE_GROUPS[1]).sum()),
                  "of_them_another_length": len(longer), "another_length_of_1_base": int((d == 1).sum()), "another_length_of_2_to_30_bases": int(((d > 1) & (d <= 30)).sum()),
                  "another_length_of_more_than_30_bases": int((d > 30).sum()), "of_them_not_determined": int((zero["change_group"] == "not determined").sum()),
                  "size_is_not_0": int((g["size"] != 0).sum())})
changed = pd.DataFrame(res)
changed.to_csv(out / "changed_bases.tsv", sep="\t", index=False)
pd.DataFrame(detail).to_csv(out / "changed_bases_detail.tsv", sep="\t", index=False)
pd.DataFrame(sizes).to_csv(out / "size_column.tsv", sep="\t", index=False)

# one substitution per event: its placements lie on the two haplotypes and are the same change
sub = S[S["change_group"] == "one changed base"]
per_event = sub.groupby(["run", "class", "event"], sort=False).agg(
    substitution=("substitution", lambda v: v.iloc[0] if v.nunique() == 1 else "differs"),
    base_before=("base_before", lambda v: v.iloc[0] if v.nunique() == 1 else "differs"),
    base_after=("base_after", lambda v: v.iloc[0] if v.nunique() == 1 else "differs")).reset_index()
res, nb = [], []
for (key, cls), g in sub.groupby(["run", "class"], sort=False):
    e = per_event[(per_event["run"] == key) & (per_event["class"] == cls)]
    e_ok = e[e["substitution"] != "differs"]
    for name in KINDS:
        x, ex = g[g["substitution"] == name], e_ok[e_ok["substitution"] == name]
        row = {"run": key, "class": cls, "substitution": name, "placements": len(x), "of": len(g), "events": len(ex), "of_events": len(e_ok),
               "share_of_events": share(len(ex), len(e_ok))}
        if name == "C>T":
            row["events_before_a_g"] = int((ex["base_after"] == "G").sum())
        res.append(row)
        for before in "ACGT":
            for after in "ACGT":
                nb.append({"run": key, "class": cls, "substitution": name, "base_before": before, "base_after": after,
                           "placements": int(((x["base_before"] == before) & (x["base_after"] == after)).sum()),
                           "events": int(((ex["base_before"] == before) & (ex["base_after"] == after)).sum())})
    ti = int(e_ok["substitution"].isin(["C>T", "T>C"]).sum())
    res.append({"run": key, "class": cls, "substitution": "transitions (C>T, T>C)", "placements": int(g["substitution"].isin(["C>T", "T>C"]).sum()), "of": len(g),
                "events": ti, "of_events": len(e_ok), "share_of_events": share(ti, len(e_ok))})
    res.append({"run": key, "class": cls, "substitution": "events whose placements differ", "placements": "", "of": "", "events": len(e) - len(e_ok),
                "of_events": len(e), "share_of_events": share(len(e) - len(e_ok), len(e))})
subs = pd.DataFrame(res)
subs.to_csv(out / "substitutions.tsv", sep="\t", index=False)
pd.DataFrame(nb).to_csv(out / "substitutions_with_neighbours.tsv", sep="\t", index=False)

res = []
lone = L[L["change_group"] == "one changed base"]
for (key, cls), g in lone.groupby(["run", "class"], sort=False):
    for name in P.GROUPS:
        x = g[g["event_group"] == name]
        e = x.drop_duplicates("event")
        row = {"run": key, "class": cls, "event_sorted_as": name, "events": len(e)}
        for kind in KINDS:
            row[kind] = int((e["substitution"] == kind).sum())
        row["C>T before a G"] = int(((e["substitution"] == "C>T") & (e["base_after"] == "G")).sum())
        res.append(row)
pd.DataFrame(res).to_csv(out / "substitutions_by_long_reads.tsv", sep="\t", index=False)

# the six kinds beside what the cluster table already says of each placement: read count, strands, and
# whether another assembly of the person holds the allele (the truth genome for hprc_year2, one or both
# HGSVC3 assemblies for NA19240; no other assembly rode along in the other two runs)
res = []
for (key, cls), g in sub.groupby(["run", "class"], sort=False):
    for name in KINDS + ["all six"]:
        x = g if name == "all six" else g[g["substitution"] == name]
        e = x.groupby("event")
        row = {"run": key, "class": cls, "substitution": name, "placements": len(x), "events": int(x["event"].nunique()),
               "read_count_median": float(x["read_count"].median()) if len(x) else "",
               "share_of_reads_on_the_rarer_strand_median": round(float((x["reads_on_the_rarer_strand"] / x["reads"]).median()), 2) if len(x) else ""}
        if g["held_elsewhere"].notna().any():
            held = int((e["held_elsewhere"].agg(lambda v: (v == "yes").any())).sum()) if len(x) else 0
            row.update({"events_held_by_another_assembly": held, "share_held": share(held, int(x["event"].nunique()))})
        res.append(row)
pd.DataFrame(res).to_csv(out / "substitutions_beside_the_cluster_table.tsv", sep="\t", index=False)

# ---- to the terminal ----------------------------------------------------------------------------------------------------
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)
print(pd.DataFrame(checks).T.to_string())
print(by_parent.pivot(index=["run", "class"], columns="long_reads", values="placements")[P.GROUPS].to_string())
print(events.pivot(index=["run", "class"], columns="long_reads", values="events")[P.GROUPS].to_string())
print(changed.pivot(index=["run", "class"], columns="change", values="placements")[CHANGE_GROUPS].to_string())
print(subs[subs["substitution"].isin(KINDS)].pivot(index=["run", "class"], columns="substitution", values="events")[KINDS].to_string())
