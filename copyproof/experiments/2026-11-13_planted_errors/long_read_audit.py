#!/usr/bin/env python3
"""Rung 2 with long reads: what does the long-read audit say at planted false adjacencies?

There is no long-read alignment to the planted genome, so one is made from real reads: the HiFi
reads that the public alignment places near each edit of the unedited genome are fetched and
aligned again to the planted genome (`scripts/run_planted_hifi.sh`). Four steps live here.

  long_read_audit.py regions --edits plant1.edits.tsv --fai GENOME.fai --out REGIONS.bed [--flank 50000]

      The stretches of the unedited genome whose reads are fetched: every edit and the stretch
      homologous to it on the other haplotype, each with `flank` bases on both sides.

  long_read_audit.py reads --regions REGIONS.bed --fetched FETCHED.bam --aligned ALIGNED.bam --out TABLE.tsv

      What was fetched and what the new alignment made of it: records by type in both files, reads
      that could not be aligned again (only a secondary or supplementary record lies in a fetched
      stretch), and reads whose primary alignment is now on another sequence.

  long_read_audit.py records --edits plant1.edits.tsv --expected plant1.expected.bed --bam ALIGNED.bam --out PREFIX

      PREFIX.records.tsv          per planted false adjacency, the alignment records there by type:
                                  how many run through it (and how many of those hold an insertion
                                  or a deletion of 500 bases or more at it), how many are clipped at it
      PREFIX.records_by_size.tsv  sums and medians of those counts per kind of adjacency and size of the edit

  long_read_audit.py score --edits plant1.edits.tsv --expected plant1.expected.bed --regions REGIONS.bed
                           --run PLANTED_RUN_WITH_LONG_READS [--short-run PLANTED_RUN]
                           [--truth-run TRUTH_RUN_WITH_LONG_READS] --out PREFIX

      PREFIX.adjacencies.tsv   one row per planted false adjacency: what short reads show, what the
                               audit measured, and the issue that results
      PREFIX.summary.tsv       the same, counted per kind of adjacency and stratum
      PREFIX.by_size.tsv       the same, counted per kind of adjacency and size of the edit
      PREFIX.elsewhere.tsv     every audited place inside the fetched stretches and away from the
                               edits, beside the verdict of the truth genome's run at the same place
      PREFIX.audit_issues.tsv  every issue the audit made or changed (channel S1+S4 or S4), with its
                               place: at an edit, fetched and away from the edits, or not fetched
      PREFIX.issue_changes.tsv issues by class and confidence, without and with the long reads

A planted false adjacency is the join left by a deletion, either end of an inversion, either run
of N of an inversion rejoined through N, and the junction between the two copies of a tandem
duplication (the last is not in `plant1.expected.bed`; it is derived from the edits table).
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

JOIN_CLASSES = ["FALSE_JOIN", "BASE_ERROR", "UNSUPPORTED", "UNSUPPORTED_JOIN"]
ADJACENCY = {"del": "deletion join", "inv": "inversion breakpoint", "invgap": "join through 100 N",
             "dup": "junction of a tandem duplication"}


def read_edits(path) -> pd.DataFrame:
    ed = pd.read_csv(path, sep="\t")
    ed["delta"] = (ed["new_end"] - ed["new_start"]) - (ed["end"] - ed["start"])
    ed.loc[ed["kind"] == "extra", "delta"] = 0  # the copy is a sequence of its own; nothing shifts
    return ed


class Coordinates:
    """Between the unedited genome and the planted one, for positions outside every edit."""

    def __init__(self, ed: pd.DataFrame):
        self.new, self.old = {}, {}
        for c, g in ed.groupby("contig"):
            g = g.sort_values("start")
            cum = np.cumsum(g["delta"].values)
            self.new[c] = (g["new_end"].values, cum)
            self.old[c] = (g["end"].values, cum)

    def to_original(self, c: str, pos: int) -> int:
        t = self.new.get(c)
        if t is None:
            return int(pos)
        i = int(np.searchsorted(t[0], pos, side="right"))
        return int(pos - (t[1][i - 1] if i else 0))

    def to_planted(self, c: str, pos: int) -> int:
        t = self.old.get(c)
        if t is None:
            return int(pos)
        i = int(np.searchsorted(t[0], pos, side="right"))
        return int(pos + (t[1][i - 1] if i else 0))


def merge(rows):
    """Merge overlapping (contig, start, end, label) rows; labels are joined."""
    out = []
    for c, s, e, lab in sorted(rows):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
            out[-1][3].append(lab)
        else:
            out.append([c, s, e, [lab]])
    return [(c, s, e, ",".join(lab)) for c, s, e, lab in out]


def cmd_regions(a) -> int:
    ed = read_edits(a.edits)
    length = {}
    if a.fai:
        length = {t[0]: int(t[1]) for t in (line.split("\t") for line in open(a.fai) if line.strip())}
    rows = []
    for r in ed.itertuples(index=False):
        for c, s, e, what in ((r.contig, r.start, r.end, "edit"), (r.hom_contig, r.hom_start, r.hom_end, "homologue")):
            if isinstance(c, str):
                rows.append((c, max(int(s) - a.flank, 0), min(int(e) + a.flank, length.get(c, 1 << 62)), f"{r.id}:{what}"))
    merged = merge(rows)
    with open(a.out, "w") as f:
        for c, s, e, lab in merged:
            f.write(f"{c}\t{s}\t{e}\t{lab}\n")
    total = sum(e - s for _, s, e, _ in merged)
    print(f"{len(ed)} edits, {len(rows)} stretches of the unedited genome, {len(merged)} after merging, {total} bases; "
          f"flank {a.flank}", file=sys.stderr)
    return 0


_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def records(bam: str):
    """(name, flag, sequence, start, end, mapping quality, read length by the CIGAR) of every record of a BAM."""
    proc = subprocess.Popen(f"samtools view '{bam}' | cut -f 1-6", shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        name, flag, contig, pos, mapq, cigar = line.rstrip("\n").split("\t")
        ops = _CIG.findall(cigar)
        start = int(pos) - 1
        yield (name, int(flag), contig, start, start + sum(int(n) for n, o in ops if o in "MDN=X"), int(mapq),
               sum(int(n) for n, o in ops if o in "MIS=XH"))
    if proc.wait():
        raise RuntimeError(f"samtools view failed on {bam}")


def cmd_reads(a) -> int:
    reg = pd.read_csv(a.regions, sep="\t", names=["contig", "start", "end", "label"])
    rows = [("fetched stretches of the unedited genome", len(reg)), ("... bases", int((reg["end"] - reg["start"]).sum()))]
    n = {"primary": 0, "secondary": 0, "supplementary": 0}
    home, other, bases, inside, q20 = {}, set(), 0, 0, 0
    rby = {c: g[["start", "end"]].values for c, g in reg.groupby("contig")}
    for name, flag, contig, start, end, mapq, length in records(a.fetched):
        kind = "secondary" if flag & 0x100 else ("supplementary" if flag & 0x800 else "primary")
        n[kind] += 1
        if kind == "primary":
            home[name] = contig
            bases += length
            q20 += mapq >= 20
            r = rby.get(contig)  # the part of the alignment that lies inside a fetched stretch
            if r is not None:
                inside += int(np.maximum(np.minimum(r[:, 1], end) - np.maximum(r[:, 0], start), 0).sum())
        else:
            other.add(name)
    rows += [("public alignment: records in those stretches", sum(n.values())),
             ("... primary (the reads that are aligned again)", n["primary"]),
             ("... ... with mapping quality 20 or more", q20),
             ("... ... bases", bases), ("... ... mean read length", round(bases / max(n["primary"], 1))),
             ("... ... mean depth of their alignments inside the fetched stretches",
              round(inside / max(int((reg["end"] - reg["start"]).sum()), 1), 1)),
             ("... secondary", n["secondary"]), ("... supplementary", n["supplementary"]),
             ("... reads with only a secondary or supplementary record there (not aligned again)", len(other - set(home)))]
    m = {"primary": 0, "secondary": 0, "supplementary": 0, "unmapped": 0}
    moved, q20, seen = 0, 0, set()
    for name, flag, contig, _, _, mapq, length in records(a.aligned):
        kind = "unmapped" if flag & 0x4 else ("secondary" if flag & 0x100 else ("supplementary" if flag & 0x800 else "primary"))
        m[kind] += 1
        if kind == "primary":
            seen.add(name)
            q20 += mapq >= 20
            moved += home.get(name) != contig
    rows += [("new alignment to the planted genome: records", sum(m.values())),
             ("... primary", m["primary"]), ("... ... with mapping quality 20 or more", q20),
             ("... ... on another sequence than in the public alignment", moved),
             ("... secondary", m["secondary"]), ("... supplementary", m["supplementary"]), ("... unmapped reads", m["unmapped"]),
             ("... reads of the public alignment that are missing from the new one", len(set(home) - seen) - m["unmapped"])]
    out = pd.DataFrame([(what, f"{n:g}" if isinstance(n, float) else str(n)) for what, n in rows], columns=["what", "n"])
    out.to_csv(a.out, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


def num(x) -> int:
    """A count from a table cell; an empty cell (a place the alignment file could not answer) is 0."""
    return 0 if pd.isna(x) else int(x)


def touching(table: pd.DataFrame, by: dict, c: str, lo: int, hi: int) -> pd.DataFrame:
    """Rows of one sequence that reach into [lo, hi], both ends included (as scripts/eval_planted.py compares)."""
    g = by.get(c)
    if g is None:
        return table.iloc[:0]
    return g[(g["start"] <= hi) & (g["end"] >= lo)]


def read_expected(path) -> pd.DataFrame:
    exp = pd.read_csv(path, sep="\t", names=["contig", "start", "end", "signal", "tag"])
    exp[["id", "kind", "size", "stratum"]] = exp["tag"].str.split("|", expand=True)
    return exp


def adjacencies(exp: pd.DataFrame, ed: pd.DataFrame) -> list:
    """The planted false adjacencies, in planted coordinates: (id, kind, stratum, size, contig, start, end).
    A join is the two bases around it; a join through N is the run of N."""
    adj = []
    for r in exp[exp["signal"].isin(["FALSE_JOIN", "MISSING_JUNCTION"])].itertuples(index=False):
        adj.append((r.id, r.kind, r.stratum, int(r.size), r.contig, int(r.start), int(r.end)))
    for r in ed[ed["kind"] == "dup"].itertuples(index=False):
        mid = int(r.new_start) + int(r.end - r.start)  # first base of the second copy
        adj.append((r.id, "dup", r.stratum, int(r.size), r.contig, mid - 1, mid + 1))
    return adj


def cmd_records(a) -> int:
    """The alignment records at every planted false adjacency, by type. The counts are taken at
    the adjacency itself with the audit's distances (joins.py, audit_points): an alignment runs
    through when it covers 500 bases on both sides, and is clipped there when it stops within
    100 bases with a clip of 200 or more."""
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    ed = read_edits(a.edits)
    homologue = dict(zip(ed["id"], ed["hom_contig"]))  # the sequence of the other haplotype, per edit
    home = {r[0]: r[2] for r in records(a.bam) if not r[1] & 0x904}  # read -> sequence of its primary record
    rows = []
    for pid, kind, stratum, size, c, s, e in adjacencies(read_expected(a.expected), ed):
        out = subprocess.run(["samtools", "view", "-F", "0x4", a.bam, f"{c}:{max(s - 2000, 0) + 1}-{e + 2000}"],
                             capture_output=True, text=True, check=True).stdout
        n = {k: 0 for k in ("records", "primary", "secondary", "supplementary", "primary_over", "through", "through_primary",
                            "through_secondary", "through_secondary_primary_on_homologue", "through_mapq20",
                            "through_large_insertion", "through_large_deletion", "clipped", "clipped_primary",
                            "clipped_secondary", "clipped_secondary_primary_on_homologue", "clipped_supplementary")}
        indels = []
        for line in out.splitlines():
            t = line.split("\t", 6)
            elsewhere = home.get(t[0]) == homologue[pid]  # this read's primary record is on the other haplotype
            flag, pos, mapq = int(t[1]), int(t[3]) - 1, int(t[4])
            ops = [(int(x), o) for x, o in _CIG.findall(t[5])]
            lead = ops[0][0] if ops and ops[0][1] in "SH" else 0
            trail = ops[-1][0] if len(ops) > 1 and ops[-1][1] in "SH" else 0
            ref, big, what = pos, 0, ""  # the largest insertion or deletion within 100 bases of the adjacency
            for x, o in ops:
                if o in "ID" and x >= a.min_indel and x > big and ref <= e + 100 and ref + (x if o == "D" else 0) >= s - 100:
                    big, what = x, "insertion" if o == "I" else "deletion"
                if o in "MDN=X":
                    ref += x
            typ = "secondary" if flag & 0x100 else ("supplementary" if flag & 0x800 else "primary")
            n["records"] += 1
            n[typ] += 1
            n["primary_over"] += typ == "primary" and pos < e and ref > s
            if pos <= s - 500 and ref >= e + 500:
                n["through"] += 1
                n["through_mapq20"] += mapq >= 20
                if typ != "supplementary":
                    n[f"through_{typ}"] += 1
                n["through_secondary_primary_on_homologue"] += typ == "secondary" and elsewhere
                if big:
                    n[f"through_large_{what}"] += 1
                    indels.append(big)
            if ((abs(ref - s) <= 100 or abs(ref - e) <= 100) and trail >= 200) or \
                    ((abs(pos - s) <= 100 or abs(pos - e) <= 100) and lead >= 200):
                n["clipped"] += 1
                n[f"clipped_{typ}"] += 1
                n["clipped_secondary_primary_on_homologue"] += typ == "secondary" and elsewhere
        rows.append({"id": pid, "kind": kind, "adjacency": ADJACENCY[kind], "stratum": stratum, "size": size, "contig": c,
                     "start": s, "end": e, **n,
                     "median_length_of_those": int(np.median(indels)) if indels else ""})
    res = pd.DataFrame(rows)
    res.to_csv(f"{a.out}.records.tsv", sep="\t", index=False)
    cols = ["primary_over", "through", "through_primary", "through_secondary", "through_secondary_primary_on_homologue",
            "through_mapq20", "through_large_insertion", "through_large_deletion", "clipped", "clipped_primary",
            "clipped_secondary", "clipped_secondary_primary_on_homologue", "clipped_supplementary"]
    # per kind and size: the number of adjacencies, then the sum of every count over them, then the medians
    g = res.groupby(["adjacency", "size"], sort=False)
    tab = pd.concat([g.size().rename("planted"), g[cols].sum().add_prefix("sum_"), g[cols].median().add_prefix("median_")],
                    axis=1).reset_index()
    tab.to_csv(f"{a.out}.records_by_size.tsv", sep="\t", index=False)
    print("alignment records at the planted false adjacencies, per kind and size (sums over the adjacencies, then medians):")
    print(tab.T.to_string(header=False))
    return 0


def cmd_score(a) -> int:
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    pd.set_option("display.max_columns", 50)
    ed = read_edits(a.edits)
    xy = Coordinates(ed)
    exp = read_expected(a.expected)
    run = Path(a.run)
    iss = pd.read_csv(run / "issues.tsv", sep="\t")
    pts = pd.read_csv(run / "unsupported_points.tsv", sep="\t", dtype={"name": str})
    unsup = pd.read_csv(run / "unsupported_clusters.tsv", sep="\t")
    joins = pd.read_csv(run / "joins.tsv", sep="\t") if (run / "joins.tsv").is_file() else pd.DataFrame()
    runs = pts[pts["name"].str.startswith("run")]
    iby = {c: g for c, g in iss.groupby("contig")}
    pby = {c: g for c, g in runs.groupby("contig")}
    uby = {c: g for c, g in unsup.groupby("contig")}

    # 1. the planted false adjacencies
    adj = adjacencies(exp, ed)

    rows = []
    tol = a.point_tol
    for pid, kind, stratum, size, c, s, e in adj:
        row = {"id": pid, "kind": kind, "adjacency": ADJACENCY[kind], "stratum": stratum, "size": size, "contig": c,
               "start": s, "end": e}
        mid = (s + e) // 2
        lo, hi = mid - tol, mid + tol
        # short reads: a run of assembly k-mers that no read holds (as scripts/eval_planted.py counts them)
        u = touching(unsup, uby, c, lo, hi)
        row["unsupported_run"] = len(u) > 0
        row["covered_flanks"] = bool(len(u) and (u["flank_copies"] >= 0.85).any())
        if len(u):  # the numbers that decide whether the pipeline audits the run (pipeline.py, step 9)
            b = u.sort_values("unsupported", ascending=False).iloc[0]
            row.update({"run_unsupported": int(b["unsupported"]), "run_tracked": int(b["tracked"]),
                        "run_flank_copies": float(b["flank_copies"]), "run_flank_kmers": int(b["flank_kmers"])})
        # the audit of that run
        p = touching(runs, pby, c, lo, hi)
        row["audited"] = len(p) > 0
        row["by_alignment"] = row["by_content"] = ""
        if len(p):
            p = p.sort_values("kmers", ascending=False).iloc[0]
            row["point"] = p["name"]
            for col in ("reads", "long_reads_over", "spanning", "spanning_mapq", "clipped", "kmers", "kmers_carried"):
                row[col] = num(p[col])
            row["by_alignment"], row["by_content"] = p["verdict"], p["long_read_kmers"]
        row["join_verdict"] = ""
        # the audit of a run of N
        if kind == "invgap" and len(joins):
            j = joins[(joins["contig"] == c) & (joins["gap_start"] < e) & (joins["gap_end"] > s)]
            if len(j):
                j = j.iloc[0]
                row["join_verdict"] = j["verdict"]
                for col in ("linking_reads", "depth_left", "depth_right", "clips_left", "clips_right"):
                    row[f"join_{col}"] = int(j[col])
        # the issues that result
        g = touching(iss, iby, c, lo, hi)
        g = g[g["class"].isin(JOIN_CLASSES)]
        order = {k: i for i, k in enumerate(JOIN_CLASSES)}
        g = g.assign(_o=g["class"].map(order)).sort_values(["_o", "confidence"])
        row["issue"] = g.iloc[0]["class"] if len(g) else "none"
        row["confidence"] = g.iloc[0]["confidence"] if len(g) else ""
        row["channel"] = g.iloc[0]["channel"] if len(g) else ""
        row["all_issues"] = ",".join(f"{k}/{f}/{ch}" for k, f, ch in zip(g["class"], g["confidence"], g["channel"]))
        rows.append(row)
    first = ["id", "kind", "adjacency", "stratum", "size", "contig", "start", "end", "unsupported_run", "covered_flanks", "run_unsupported",
             "run_tracked", "run_flank_copies", "run_flank_kmers", "audited", "point", "reads", "long_reads_over", "spanning", "spanning_mapq", "clipped", "kmers", "kmers_carried", "by_alignment",
             "by_content", "join_verdict", "join_linking_reads", "join_depth_left", "join_depth_right", "join_clips_left",
             "join_clips_right", "issue", "confidence", "channel", "all_issues"]
    res = pd.DataFrame(rows).reindex(columns=first)
    for col in ("run_unsupported", "run_tracked", "run_flank_kmers", "reads", "long_reads_over", "spanning", "spanning_mapq",
                "clipped", "kmers", "kmers_carried", "join_linking_reads", "join_depth_left", "join_depth_right",
                "join_clips_left", "join_clips_right"):
        if col in res:
            res[col] = res[col].astype("Int64")
    res.to_csv(f"{a.out}.adjacencies.tsv", sep="\t", index=False)

    def count(g, cls, conf=None):
        m = g["issue"] == cls
        if conf:
            m &= g["confidence"] == conf
        return int(m.sum())

    def median(col):
        return "" if col.dropna().empty else f"{float(col.dropna().median()):g}"

    def tally(g, **key):
        au = g[g["audited"]]
        return {
            **key, "planted": len(g),
            "unsupported_run": int(g["unsupported_run"].sum()),
            "run_with_covered_flanks": int(g["covered_flanks"].sum()), "audited": int(g["audited"].sum()),
            "median_primary_reads_over": median(au["long_reads_over"]), "median_spanning": median(au["spanning"]),
            "median_clipped": median(au["clipped"]),
            "some_unread_kmer_in_a_long_read": int((au["kmers_carried"] > 0).sum()),
            "content_not_carried": int((g["by_content"] == "not_carried").sum()),
            "content_carried": int((g["by_content"] == "carried").sum()),
            "content_unclear": int((g["by_content"] == "unclear").sum()),
            "alignment_broken": int((g["by_alignment"] == "broken").sum()),
            "alignment_spanned": int((g["by_alignment"] == "spanned").sum()),
            "alignment_unclear": int((g["by_alignment"] == "unclear").sum()),
            "join_unsupported": int((g["join_verdict"] == "unsupported").sum()),
            "join_linked": int((g["join_verdict"] == "linked").sum()),
            "join_not_assessable": int((g["join_verdict"] == "not_assessable").sum()),
            "FALSE_JOIN_high": count(g, "FALSE_JOIN", "high"), "FALSE_JOIN_low": count(g, "FALSE_JOIN", "low"),
            "BASE_ERROR_high": count(g, "BASE_ERROR", "high"), "BASE_ERROR_low": count(g, "BASE_ERROR", "low"),
            "UNSUPPORTED": count(g, "UNSUPPORTED"),
            "UNSUPPORTED_JOIN_high": count(g, "UNSUPPORTED_JOIN", "high"),
            "UNSUPPORTED_JOIN_low": count(g, "UNSUPPORTED_JOIN", "low"),
            "no_issue": count(g, "none"),
        }

    summ, by_size = [], []
    for kind in ("inv", "del", "invgap", "dup"):
        k = res[res["kind"] == kind]
        strata = sorted(k["stratum"].unique(), reverse=True)
        for stratum in (strata + ["all"] if len(strata) > 1 else ["all"]):
            summ.append(tally(k[(k["stratum"] == stratum) | (stratum == "all")], adjacency=ADJACENCY[kind], stratum=stratum))
        for size in sorted(k["size"].unique()):
            by_size.append(tally(k[k["size"] == size], adjacency=ADJACENCY[kind], size=int(size)))
    summ = pd.DataFrame(summ)
    summ.to_csv(f"{a.out}.summary.tsv", sep="\t", index=False)
    pd.DataFrame(by_size).to_csv(f"{a.out}.by_size.tsv", sep="\t", index=False)
    print("planted false adjacencies, by what the long-read audit made of them:")
    print(summ.T.to_string(header=False))
    print("\nthe same by size of the planted edit:")
    print(pd.DataFrame(by_size).T.to_string(header=False))

    # 2. inside the fetched stretches, away from every edit
    reg = pd.read_csv(a.regions, sep="\t", names=["contig", "start", "end", "label"])
    rby = {c: g for c, g in reg.groupby("contig")}
    eby = {c: g[["new_start", "new_end"]].values for c, g in ed.groupby("contig")}
    for r in ed[ed["kind"] == "extra"].itertuples(index=False):  # the extra sequences are edits as a whole
        eby[f"{r.id}_extra"] = np.array([[0, int(r.end - r.start)]])
    hby = {c: g[["hom_start", "hom_end"]].values for c, g in ed.dropna(subset=["hom_contig"]).groupby("hom_contig")}

    def edit_distance(c, s, e):
        """Bases between a stretch (planted coordinates) and the nearest edit of its sequence."""
        t = eby.get(c)
        if t is None:
            return np.inf
        return float(np.min(np.maximum(np.maximum(t[:, 0] - e, s - t[:, 1]), 0)))

    def fetched(c, s, e):
        """The fetched stretch (original coordinates) that holds this stretch of the planted genome."""
        g = rby.get(c)
        if g is None:
            return ""
        so, eo = xy.to_original(c, s), xy.to_original(c, e)
        g = g[(g["start"] <= so) & (g["end"] >= eo)]
        return g.iloc[0]["label"] if len(g) else ""

    def homologue_distance(c, s, e):
        t = hby.get(c)
        if t is None:
            return np.inf
        so, eo = xy.to_original(c, s), xy.to_original(c, e)
        return float(np.min(np.maximum(np.maximum(t[:, 0] - eo, so - t[:, 1]), 0)))

    truth = None
    if a.truth_run:
        truth = pd.read_csv(Path(a.truth_run) / "unsupported_points.tsv", sep="\t", dtype={"name": str})
        truth = {c: g for c, g in truth[truth["name"].str.startswith("run")].groupby("contig")}
    else_rows = []
    for p in runs.itertuples(index=False):
        d = edit_distance(p.contig, p.start, p.end)
        lab = fetched(p.contig, p.start, p.end)
        if not lab or d <= a.pad:
            continue
        row = {"name": p.name, "contig": p.contig, "start": int(p.start), "end": int(p.end),
               "original_start": xy.to_original(p.contig, int(p.start)), "fetched_for": lab,
               "distance_to_edit": int(d) if np.isfinite(d) else -1,
               "distance_to_homologue": (lambda h: int(h) if np.isfinite(h) else -1)(homologue_distance(p.contig, p.start, p.end)),
               "reads": num(p.reads), "long_reads_over": num(p.long_reads_over), "spanning": num(p.spanning),
               "clipped": num(p.clipped), "kmers": num(p.kmers), "kmers_carried": num(p.kmers_carried),
               "by_alignment": p.verdict, "by_content": p.long_read_kmers}
        g = touching(iss, iby, p.contig, int(p.start), int(p.end))
        g = g[g["class"].isin(JOIN_CLASSES) & (g["channel"] == "S1+S4")]
        row["issue"] = g.iloc[0]["class"] if len(g) else "none"
        row["confidence"] = g.iloc[0]["confidence"] if len(g) else ""
        if truth is not None:
            so = row["original_start"]
            t = touching(runs, truth, p.contig, so, so + int(p.end - p.start))
            if len(t):
                t = t.iloc[0]
                row.update({"truth_reads": num(t["reads"]), "truth_long_reads_over": num(t["long_reads_over"]),
                            "truth_spanning": num(t["spanning"]), "truth_clipped": num(t["clipped"]), "truth_kmers": num(t["kmers"]),
                            "truth_kmers_carried": num(t["kmers_carried"]), "truth_by_alignment": t["verdict"],
                            "truth_by_content": t["long_read_kmers"]})
            else:
                row.update({"truth_by_alignment": "no_such_run", "truth_by_content": "no_such_run"})
        else_rows.append(row)
    away = pd.DataFrame(else_rows)
    for col in [c for c in away.columns if c.startswith("truth_") and not c.startswith("truth_by")]:
        away[col] = away[col].astype("Int64")
    away.to_csv(f"{a.out}.elsewhere.tsv", sep="\t", index=False)
    in_reg = [bool(fetched(c, s, e)) for c, s, e in zip(iss["contig"], iss["start"], iss["end"])]
    dist = np.array([edit_distance(c, s, e) for c, s, e in zip(iss["contig"], iss["start"], iss["end"])])
    iss["where"] = np.where(dist <= a.pad, "at an edit", np.where(in_reg, "fetched, away from the edits", "not fetched"))
    print(f"\naudited runs inside the fetched stretches and more than {a.pad} bp from every edit: {len(away)}")
    if len(away):
        print(away.groupby(["by_content", "by_alignment", "issue"]).size().to_string())
        if truth is not None:
            print("\nthe same runs in the truth genome's run (public alignment), content verdict there against here:")
            print(pd.crosstab(away["truth_by_content"], away["by_content"]).to_string())
            print("alignment verdict there against here:")
            print(pd.crosstab(away["truth_by_alignment"], away["by_alignment"]).to_string())
    iss["distance_to_edit"] = np.where(np.isfinite(dist), dist, -1).astype(int)
    made = iss[iss["channel"].isin(["S1+S4", "S4"])]  # what the audit made or changed
    made.drop(columns=["summary"]).to_csv(f"{a.out}.audit_issues.tsv", sep="\t", index=False)
    print("\nissues that the long-read audit made or changed, by place:")
    print(made.groupby(["where", "class", "channel", "confidence"]).size().to_string())
    s4 = iss[iss["class"].isin(JOIN_CLASSES)]
    print("\nall issues of the join classes, by place:")
    print(s4.groupby(["where", "class", "channel", "confidence"]).size().to_string())

    # 3. what the long reads changed in the issue list
    if a.short_run:
        short = pd.read_csv(Path(a.short_run) / "issues.tsv", sep="\t")
        t = pd.concat([short.groupby(["class", "confidence"]).size().rename("short_reads_only"),
                       iss.groupby(["class", "confidence"]).size().rename("with_long_reads")], axis=1).fillna(0).astype(int)
        t["difference"] = t["with_long_reads"] - t["short_reads_only"]
        t.reset_index().to_csv(f"{a.out}.issue_changes.tsv", sep="\t", index=False)
        print(f"\nissues: {len(short)} with short reads only, {len(iss)} with long reads")
        print(t[t["difference"] != 0].to_string())
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("regions")
    r.add_argument("--edits", required=True)
    r.add_argument("--fai", help="index of the unedited genome (sequence lengths)")
    r.add_argument("--flank", type=int, default=50_000)
    r.add_argument("--out", required=True)
    q = sub.add_parser("reads")
    q.add_argument("--regions", required=True)
    q.add_argument("--fetched", required=True, help="the records of the public alignment in the fetched stretches")
    q.add_argument("--aligned", required=True, help="the same reads aligned to the planted genome")
    q.add_argument("--out", required=True)
    c = sub.add_parser("records")
    c.add_argument("--edits", required=True)
    c.add_argument("--expected", required=True)
    c.add_argument("--bam", required=True, help="the long reads aligned to the planted genome")
    c.add_argument("--min-indel", type=int, default=500, help="an insertion or deletion this long counts as large")
    c.add_argument("--out", required=True)
    s = sub.add_parser("score")
    s.add_argument("--edits", required=True)
    s.add_argument("--expected", required=True)
    s.add_argument("--regions", required=True)
    s.add_argument("--run", required=True)
    s.add_argument("--short-run")
    s.add_argument("--truth-run")
    s.add_argument("--out", required=True)
    s.add_argument("--point-tol", type=int, default=200, help="as in scripts/eval_planted.py")
    s.add_argument("--pad", type=int, default=5000, help="as in scripts/eval_planted.py: nearer than this is at the edit")
    a = ap.parse_args()
    return {"regions": cmd_regions, "reads": cmd_reads, "records": cmd_records, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
