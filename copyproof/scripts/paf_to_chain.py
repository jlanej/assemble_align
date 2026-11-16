#!/usr/bin/env python3
"""Turn an alignment of one haplotype against the other into a chain file for scripts/plant_errors.py.

The planting tool needs, for every place it edits, the homologous stretch of the other haplotype.
A T2T genome ships a chain between its haplotypes; for any other diploid release, align the two
haplotype FASTAs and convert:

  minimap2 -cx asm5 -t 10 pat.fa.gz mat.fa.gz > pat_mat.paf        # PAF with CIGAR (-c)
  paf_to_chain.py pat_mat.paf > hap_to_hap.chain

Only primary alignments on the plus strand are used (the planting tool lifts through plus-strand
chains only), at least --min-length long and of mapping quality at least --min-mapq. Alignments
are taken longest first, and one that overlaps an alignment already taken, on either sequence, is
left out: every position then lifts to at most one place. A position in an alignment that was left
out does not lift, and nothing is planted there. Each alignment is written twice, once from each
haplotype to the other.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def blocks(cigar: str):
    """(size, gap on the target, gap on the query) per aligned block; the last block has no gaps."""
    out, size, dt, dq = [], 0, 0, 0
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            if dt or dq:
                out.append((size, dt, dq))
                size, dt, dq = 0, 0, 0
            size += n
        elif op == "I":
            dq += n
        elif op in "DN":
            dt += n
    out.append((size, 0, 0))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf")
    ap.add_argument("--min-length", type=int, default=50_000, help="shortest alignment used, in target bases")
    ap.add_argument("--min-mapq", type=int, default=60)
    a = ap.parse_args()

    alns = []
    seen = kept_strand = 0
    with open(a.paf) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            seen += 1
            tags = {x[:2]: x[5:] for x in t[12:]}
            if tags.get("tp", "P") != "P" or t[4] != "+" or "cg" not in tags:
                continue
            kept_strand += 1
            qs, qe, ts, te, mapq = int(t[2]), int(t[3]), int(t[7]), int(t[8]), int(t[11])
            if te - ts < a.min_length or mapq < a.min_mapq:
                continue
            alns.append((te - ts, t[5], int(t[6]), ts, te, t[0], int(t[1]), qs, qe, int(t[9]), tags["cg"]))
    alns.sort(key=lambda x: -x[0])
    taken = defaultdict(list)

    def free(name, s, e):
        return all(e <= x or s >= y for x, y in taken[name])

    n = bases = 0
    for _, tn, tl, ts, te, qn, ql, qs, qe, score, cg in alns:
        if not (free(tn, ts, te) and free(qn, qs, qe)):
            continue
        taken[tn].append((ts, te))
        taken[qn].append((qs, qe))
        bl = blocks(cg)
        assert sum(b[0] + b[1] for b in bl) == te - ts and sum(b[0] + b[2] for b in bl) == qe - qs, "CIGAR and coordinates disagree"
        n += 1
        bases += te - ts
        for flip in (False, True):
            if flip:
                head = f"chain {score} {qn} {ql} + {qs} {qe} {tn} {tl} + {ts} {te} {2 * n}"
            else:
                head = f"chain {score} {tn} {tl} + {ts} {te} {qn} {ql} + {qs} {qe} {2 * n - 1}"
            out = [head]
            for size, dt, dq in bl[:-1]:
                out.append(f"{size}\t{dq}\t{dt}" if flip else f"{size}\t{dt}\t{dq}")
            out.append(str(bl[-1][0]))
            sys.stdout.write("\n".join(out) + "\n\n")
    print(f"{seen} alignments read, {kept_strand} primary on the plus strand, {len(alns)} long enough and of mapping quality "
          f"{a.min_mapq} or more, {n} used without overlap ({bases / 1e9:.3f} Gb of the target)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
