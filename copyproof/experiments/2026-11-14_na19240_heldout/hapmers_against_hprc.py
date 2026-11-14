#!/usr/bin/env python3
"""The hapmer audit on an assembly phased without the parents, held against trio-phased haplotypes
of the same person (NOTEBOOK entry p, "a run added by the user").

The run's scan carries the exact k-mers of the two trio-phased haplotypes as label sets. A tracked
k-mer of the assembly that only one of them holds marks the parent of the stretch it lies in, by an
assembly that was phased another way. For every call of the audit (phase switch or misplaced
piece), and for every 50-kb bin of the assessed sequences, the two kinds of k-mer are counted:

  call confirmed       k-mers of the parent the call names: at least MIN, and at least RATIO times
                       those of the sequence's own parent
  call contradicted    the reverse
  bin follows the      k-mers of the other parent than the sequence's: at least MIN, and at least
  other parent         RATIO times those of the sequence's own parent
  bin flagged          calls cover at least half of it

  hapmers_against_hprc.py RUN_DIR OUT_PREFIX --labels mother=hprc_mat father=hprc_pat

Writes OUT_PREFIX.sequences.tsv, .calls.tsv, .other_parent_bins.tsv and .summary.tsv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

here = Path(__file__).resolve()
sys.path.insert(0, str(here.parents[2] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402
from copyproof.track import Scan  # noqa: E402

MIN, RATIO, BIN = 20, 3.0, 50_000


def verdict(named: int, own: int) -> str:
    if named >= MIN and named >= RATIO * own:
        return "confirmed"
    if own >= MIN and own >= RATIO * named:
        return "contradicted"
    return "unclear"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--labels", nargs=2, required=True, metavar="PARENT=LABEL",
                    help="for each parent of the run, the label set of that parent's trio-phased haplotype")
    a = ap.parse_args()
    run = Path(a.run)
    sc = Scan(str(run / "scan"))
    label_of = dict(x.split("=", 1) for x in a.labels)
    p1, p2 = sc.reads[1], sc.reads[2]
    idx = {p: sc.reads.index(label_of[p]) for p in (p1, p2)}
    other = {p1: p2, p2: p1}

    def table(name, columns):
        try:  # a release without gaps has no pieces, and the audit then writes an empty file
            return pd.read_csv(run / name, sep="\t")
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=columns)

    seqs = pd.read_csv(run / "hapmer_sequences.tsv", sep="\t")
    sw = table("hapmer_switches.tsv", ["contig", "start", "end", "scaffold_parent", "segment_parent", "sites_segment_parent",
                                       "sites_scaffold_parent", "llr"])
    pc = table("hapmer_pieces.tsv", ["contig", "piece", "start", "end", "scaffold_parent", "sites_scaffold_parent",
                                     "sites_other_parent", "call"])
    calls = [sw.assign(kind="PHASE_SWITCH")[["kind", "contig", "start", "end", "scaffold_parent", "segment_parent",
                                             "sites_segment_parent", "sites_scaffold_parent"]]]
    if len(pc):
        o = pc[pc["call"] == "other_parent"]
        calls.append(pd.DataFrame({
            "kind": "MISPLACED_PIECE", "contig": o["contig"], "start": o["start"], "end": o["end"],
            "scaffold_parent": o["scaffold_parent"], "segment_parent": o["scaffold_parent"].map(other),
            "sites_segment_parent": o["sites_other_parent"], "sites_scaffold_parent": o["sites_scaffold_parent"]}))
    calls = pd.concat(calls, ignore_index=True)
    called = IntervalSet.from_rows(zip(calls["contig"], calls["start"], calls["end"]))
    by_contig = {c: g for c, g in calls.groupby("contig")}

    seq_rows, call_rows, bin_rows = [], [], []
    bins_total = bins_other = bins_other_flagged = bins_flagged = bins_flagged_other = 0
    for s in seqs.itertuples(index=False):
        pos, m, c = sc.records(s.contig)
        only = {p: (c[:, idx[p]] > 0) & (c[:, idx[other[p]]] == 0) for p in (p1, p2)}
        n = {p: int(only[p].sum()) for p in (p1, p2)}
        by_labels = p1 if verdict(n[p1], n[p2]) == "confirmed" else (p2 if verdict(n[p2], n[p1]) == "confirmed" else "unclear")
        seq_rows.append({"contig": s.contig, "group": s.group, "length": s.length, "parent_by_hapmers": s.parent,
                         f"kmers_{p1}_haplotype": n[p1], f"kmers_{p2}_haplotype": n[p2], "parent_by_labels": by_labels})
        cum = {p: np.concatenate([[0], np.cumsum(only[p])]) for p in (p1, p2)}
        for r in by_contig.get(s.contig, calls.iloc[:0]).itertuples(index=False):
            i, j = np.searchsorted(pos, [r.start, r.end])
            named = int(cum[r.segment_parent][j] - cum[r.segment_parent][i])
            own = int(cum[r.scaffold_parent][j] - cum[r.scaffold_parent][i])
            call_rows.append({**r._asdict(), "kmers_segment_parent": named, "kmers_scaffold_parent": own, "verdict": verdict(named, own)})
        if s.parent not in (p1, p2):
            continue
        nb = -(-int(s.length) // BIN)
        b = pos // BIN
        own = np.bincount(b[only[s.parent]], minlength=nb)
        oth = np.bincount(b[only[other[s.parent]]], minlength=nb)
        is_other = (oth >= MIN) & (oth >= RATIO * own)
        flagged = called.window_fraction(s.contig, int(s.length), BIN) >= 0.5
        bins_total += nb
        bins_other += int(is_other.sum())
        bins_other_flagged += int((is_other & flagged).sum())
        bins_flagged += int(flagged.sum())
        bins_flagged_other += int((flagged & is_other).sum())
        for k in np.flatnonzero(is_other):
            bin_rows.append({"contig": s.contig, "start": int(k) * BIN, "end": min(int(k + 1) * BIN, int(s.length)),
                             "sequence_parent": s.parent, "kmers_sequence_parent": int(own[k]), "kmers_other_parent": int(oth[k]),
                             "flagged": bool(flagged[k])})
    sq = pd.DataFrame(seq_rows)
    cl = pd.DataFrame(call_rows)
    sq.to_csv(f"{a.out}.sequences.tsv", sep="\t", index=False)
    cl.to_csv(f"{a.out}.calls.tsv", sep="\t", index=False)
    pd.DataFrame(bin_rows, columns=["contig", "start", "end", "sequence_parent", "kmers_sequence_parent", "kmers_other_parent",
                                    "flagged"]).to_csv(f"{a.out}.other_parent_bins.tsv", sep="\t", index=False)

    assessed = sq[sq["parent_by_hapmers"].isin([p1, p2])]
    rows = [("sequences", len(sq)), ("sequences assessed by the audit", len(assessed)),
            ("... bases in them", int(assessed["length"].sum())),
            ("... parent agrees with the trio-phased haplotypes", int((assessed["parent_by_hapmers"] == assessed["parent_by_labels"]).sum())),
            ("... parent disagrees", int(((assessed["parent_by_hapmers"] != assessed["parent_by_labels"]) & (assessed["parent_by_labels"] != "unclear")).sum())),
            ("... unclear by the trio-phased haplotypes", int((assessed["parent_by_labels"] == "unclear").sum()))]
    if len(cl):
        cl["strong"] = cl["sites_segment_parent"] >= 10  # the audit's own line between high and low confidence
        for kind in ("PHASE_SWITCH", "MISPLACED_PIECE"):
            k = cl[cl["kind"] == kind]
            sets = [("", k)] + ([(" with 10 sites or more", k[k["strong"]]), (" with under 10 sites", k[~k["strong"]])] if kind == "PHASE_SWITCH" else [])
            for name, g in sets:
                rows.append((f"{kind} calls{name}", len(g)))
                rows.append((f"... bases{name}", int((g["end"] - g["start"]).sum())))
                for v in ("confirmed", "contradicted", "unclear"):
                    rows.append((f"... {v}{name}", int((g["verdict"] == v).sum())))
    rows += [("50-kb bins of assessed sequences", bins_total), ("bins that follow the other parent by the trio-phased haplotypes", bins_other),
             ("... of those, flagged by the audit", bins_other_flagged), ("bins flagged by the audit", bins_flagged),
             ("... of those, following the other parent by the trio-phased haplotypes", bins_flagged_other)]
    out = pd.DataFrame(rows, columns=["what", "n"])
    out.to_csv(f"{a.out}.summary.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
