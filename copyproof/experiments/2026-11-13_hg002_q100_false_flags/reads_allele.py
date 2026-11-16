#!/usr/bin/env python3
"""At runs of assembly k-mers that no read holds: what do the long reads hold instead?

  reads_allele.py RUNS.tsv ASSEMBLY.fa.gz BAM_OR_URL OUT.tsv [--all]

RUNS.tsv is unsupported_runs_hifi.tsv (from annotate_runs.py). For every run that long reads do
not carry (by default only those on no known issue and unchanged in v1.2; --all for every one),
the assembly's sequence between two anchors is compared with what the long reads hold between the
same anchors. The anchors are 25 bases of the assembly on either side of the run, outside every
unread k-mer, so reads can be matched without an aligner; the alignment file only says which
reads lie there. The FASTA must be bgzip-compressed and indexed (samtools faidx).
"""

from __future__ import annotations

import subprocess
import sys
from collections import Counter

import pandas as pd

K, ANCHOR, GAP = 31, 25, 8
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def describe(ref: str, alt: str) -> str:
    """The difference between two sequences in a few words."""
    if ref == alt:
        return "same"
    a = 0
    while a < min(len(ref), len(alt)) and ref[a] == alt[a]:
        a += 1
    b = 0
    while b < min(len(ref), len(alt)) - a and ref[-1 - b] == alt[-1 - b]:
        b += 1
    r, q = ref[a : len(ref) - b], alt[a : len(alt) - b]
    if len(r) == len(q) == 1:
        return f"substitution {r}>{q}"
    if len(r) == len(q):
        return f"{len(r)} bases differ ({r}>{q})"
    longer, kind = (q, "insertion") if len(q) > len(r) else (r, "deletion")
    unit = longer if not (r and q) else ""
    left = ref[max(a - 6, 0) : a]
    if unit and len(set(unit)) == 1 and left.endswith(unit[0]):
        return f"homopolymer {unit[0]}: reads have {abs(len(q) - len(r))} {'more' if kind == 'insertion' else 'fewer'}"
    return f"{kind} of {abs(len(q) - len(r))} in the reads ({r or '-'}>{q or '-'})"


def main() -> None:
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    runs_path, fasta, bam, out = args
    runs = pd.read_csv(runs_path, sep="\t")
    sel = runs[runs["long_read_kmers"] == "not_carried"]
    if "--all" not in sys.argv:
        sel = sel[~sel["known_issue"] & ~sel["fixed_in_v1.2"]]
    rows = []
    for r in sel.itertuples(index=False):
        a = max(int(r.start) - GAP - ANCHOR, 0)
        b = int(r.end) + GAP + ANCHOR  # `end` already includes the last unread k-mer
        fa = subprocess.run(["samtools", "faidx", fasta, f"{r.contig}:{a + 1}-{b}"], capture_output=True, text=True, check=True).stdout
        ref = "".join(fa.split("\n")[1:]).upper()
        left, right = ref[:ANCHOR], ref[-ANCHOR:]
        between_ref = ref[ANCHOR:-ANCHOR]
        sam = subprocess.run(["samtools", "view", "-F", "0x904", bam, f"{r.contig}:{max(int(r.start) - 200, 1)}-{int(r.end) + 200}"],
                             capture_output=True, text=True).stdout
        seen = Counter()
        reads = 0
        for line in sam.splitlines():
            seq = line.split("\t", 10)[9].upper()
            reads += 1
            for s in (seq, revcomp(seq)):
                i = s.find(left)
                j = s.find(right, i + ANCHOR) if i >= 0 else -1
                if i >= 0 and j >= 0 and j - i < 4 * len(ref):
                    seen[s[i + ANCHOR : j]] += 1
                    break
        total = sum(seen.values())
        top, n_top = seen.most_common(1)[0] if seen else ("", 0)
        rows.append({
            "contig": r.contig, "start": int(r.start), "end": int(r.end), "known_issue": bool(r.known_issue),
            "changed_in_v1.2": bool(getattr(r, "_14", False)) if False else bool(runs.loc[(runs["contig"] == r.contig) & (runs["start"] == r.start), "fixed_in_v1.2"].iloc[0]),
            "long_reads": reads, "reads_with_both_anchors": total, "reads_like_assembly": seen.get(between_ref, 0),
            "reads_with_commonest_version": n_top, "difference": describe(between_ref, top) if total else "no read holds both anchors",
            "assembly": between_ref, "reads": top,
        })
        print({k: rows[-1][k] for k in ("contig", "start", "long_reads", "reads_with_both_anchors", "reads_like_assembly",
                                        "reads_with_commonest_version", "difference")}, flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
