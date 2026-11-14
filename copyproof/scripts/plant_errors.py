#!/usr/bin/env python3
"""Plant assembly errors of known class, size and place in a truth genome (GUIDE section 8, rung 2).

The reads come from the true genome, so after planting every edit should be found and nothing
else should change. Edits are placed on autosomes, far apart, in two strata: unique sequence and
annotated segmental duplications.

  dup     duplicate a segment in tandem                      -> EXPANSION (both copies)
  extra   ship a copy of a segment as an extra sequence      -> EXPANSION (the extra copy)
  del     delete a segment from one haplotype                -> on the homologue: COLLAPSE and
                                                                MISSING_HAPLOTYPE; a false join here
  swap    replace a segment by the other haplotype's         -> MISSING_HAPLOTYPE (and PHASE_SWITCH
                                                                with parents)
  inv     invert a segment in place                          -> two MISSING_JUNCTIONs
  invgap  invert a segment and rejoin with 100 N             -> two MISSING_JUNCTIONs at gaps

Outputs: PREFIX.fa (one line per sequence), PREFIX.edits.tsv (every edit, in original and edited
coordinates) and PREFIX.expected.bed (where each signal is expected, in edited coordinates).
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from copyproof.intervals import IntervalSet  # noqa: E402

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def revcomp(b: bytes) -> bytes:
    return b.translate(COMP)[::-1]


def read_fasta(path: str) -> dict[str, bytes]:
    """Whole FASTA in memory, upper case."""
    cmd = ["gzip", "-dc", path] if path.endswith(".gz") else ["cat", path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 24)
    seqs: dict[str, list[bytes]] = {}
    name = None
    for line in proc.stdout:
        if line[:1] == b">":
            name = line[1:].split()[0].decode()
            seqs[name] = []
        else:
            seqs[name].append(line.rstrip(b"\n"))
    proc.wait()
    return {n: b"".join(parts).upper() for n, parts in seqs.items()}


class Chain:
    """Point liftover through the aligned blocks of every plus-strand chain of a sequence."""

    def __init__(self, path: str):
        blocks = defaultdict(list)  # source sequence -> (source start, target start, size, target sequence)
        opener = gzip.open if path.endswith(".gz") else open
        cur = None
        with opener(path, "rt") as f:
            for line in f:
                t = line.split()
                if not t:
                    continue
                if t[0] == "chain":
                    cur = {"t": t[2], "q": t[7], "ok": t[4] == "+" and t[9] == "+", "tp": int(t[5]), "qp": int(t[10])}
                elif cur is not None:
                    size = int(t[0])
                    if cur["ok"]:
                        blocks[cur["t"]].append((cur["tp"], cur["qp"], size, cur["q"]))
                    if len(t) == 3:
                        cur["tp"] += size + int(t[1])
                        cur["qp"] += size + int(t[2])
        self.best = {}
        for name, bl in blocks.items():
            bl.sort()
            self.best[name] = {"t": np.array([x[0] for x in bl]), "qs": np.array([x[1] for x in bl]),
                               "n": np.array([x[2] for x in bl]), "q": [x[3] for x in bl]}

    def lift(self, contig: str, pos: int):
        """(other sequence, position), or None if the position is not in an aligned block."""
        c = self.best.get(contig)
        if c is None:
            return None
        i = int(np.searchsorted(c["t"], pos, side="right")) - 1
        if i < 0 or pos >= c["t"][i] + c["n"][i]:
            return None
        return c["q"][i], int(c["qs"][i] + (pos - c["t"][i]))

    def lift_interval(self, contig: str, start: int, end: int, slack: int = 2000):
        """Lift both ends, nudging each inwards to the nearest aligned base."""
        a = b = None
        for d in range(0, slack, 50):
            a = a or self.lift(contig, start + d)
            b = b or self.lift(contig, end - 1 - d)
            if a and b:
                break
        if not a or not b or a[0] != b[0] or b[1] <= a[1]:
            return None
        return a[0], a[1], b[1] + 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chain", required=True, help="chain from each haplotype to the other")
    ap.add_argument("--sd", required=True, help="BED of segmental duplications")
    ap.add_argument("--avoid", nargs="*", default=[], help="BEDs to stay away from (satellites, known issues)")
    ap.add_argument("--assessable", help="BED of windows the copy-number channel cannot assess (kept out of placements)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--spacing", type=int, default=3_000_000)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    print("reading the genome", file=sys.stderr)
    seqs = read_fasta(a.fasta)
    chain = Chain(a.chain)
    sd = IntervalSet.from_bed(a.sd)
    avoid = IntervalSet()
    for p in a.avoid:
        avoid = avoid.union(IntervalSet.from_bed(p))
    if a.assessable:
        avoid = avoid.union(IntervalSet.from_bed(a.assessable))
    autosomes = [n for n in seqs if n.startswith("chr") and not n.startswith(("chrX", "chrY", "chrM", "chrEBV"))]

    plan = []
    for size in (2_000, 5_000, 10_000, 50_000, 200_000):
        for stratum in ("unique", "SD"):
            plan += [("dup", size, stratum)] * a.reps + [("del", size, stratum)] * a.reps
    for size in (10_000, 50_000):
        plan += [("extra", size, "unique")] * a.reps
    for size in (10_000, 50_000, 200_000):
        plan += [("swap", size, "unique")] * a.reps + [("inv", size, "unique")] * a.reps
    plan += [("invgap", 50_000, "unique")] * a.reps

    taken = defaultdict(list)  # contig -> [(start, end)] of edits and of their homologous stretches

    def free(contig, s, e):
        return all(e + a.spacing <= ts or s - a.spacing >= te for ts, te in taken[contig])

    sd_blocks = [(c, int(s), int(e)) for c in autosomes if c in sd.by for s, e in sd.by[c]]
    edits = []
    for kind, size, stratum in plan:
        for _ in range(4000):
            if stratum == "SD":
                cands = [b for b in sd_blocks if b[2] - b[1] >= size + 2000]
                if not cands:
                    break
                c, bs, be = cands[rng.integers(len(cands))]
                s = int(rng.integers(bs + 1000, be - size - 999))
            else:
                c = autosomes[rng.integers(len(autosomes))]
                L = len(seqs[c])
                s = int(rng.integers(5_000_000, L - 5_000_000 - size))
            e = s + size
            pad = 20_000
            if b"N" in seqs[c][s - 1000 : e + 1000]:
                continue
            if avoid.overlap(c, s - pad, e + pad) > 0:
                continue
            in_sd = sd.overlap(c, s, e)
            if stratum == "unique" and sd.overlap(c, s - pad, e + pad) > 0:
                continue
            if stratum == "SD" and in_sd < 0.95 * size:
                continue
            hom = chain.lift_interval(c, s, e)
            if hom is None or abs((hom[2] - hom[1]) - size) > 0.1 * size:
                continue
            if avoid.overlap(hom[0], hom[1] - pad, hom[2] + pad) > 0:
                continue
            if not free(c, s, e) or not free(hom[0], hom[1], hom[2]):
                continue
            taken[c].append((s, e))
            taken[hom[0]].append((hom[1], hom[2]))
            edits.append({"kind": kind, "size": size, "stratum": stratum, "contig": c, "start": s, "end": e,
                          "hom_contig": hom[0], "hom_start": hom[1], "hom_end": hom[2]})
            break
        else:
            print(f"no place found for {kind} {size} {stratum}", file=sys.stderr)
    print(f"placed {len(edits)} of {len(plan)} edits", file=sys.stderr)

    # apply, right to left on each sequence, and record how coordinates shift
    by_contig = defaultdict(list)
    for i, ed in enumerate(edits):
        ed["id"] = f"plant{i + 1:03d}"
        by_contig[ed["contig"]].append(ed)
    shifts = {}  # contig -> sorted [(original position, cumulative shift after it)]
    extras = {}
    new = dict(seqs)
    for c, eds in by_contig.items():
        s = new[c]
        for ed in sorted(eds, key=lambda x: -x["start"]):
            a0, b0 = ed["start"], ed["end"]
            seg = s[a0:b0]
            if ed["kind"] == "dup":
                s = s[:b0] + seg + s[b0:]
                ed["delta"] = b0 - a0
            elif ed["kind"] == "del":
                s = s[:a0] + s[b0:]
                ed["delta"] = -(b0 - a0)
            elif ed["kind"] == "inv":
                s = s[:a0] + revcomp(seg) + s[b0:]
                ed["delta"] = 0
            elif ed["kind"] == "invgap":
                s = s[:a0] + b"N" * 100 + revcomp(seg) + b"N" * 100 + s[b0:]
                ed["delta"] = 200
            elif ed["kind"] == "swap":
                other = seqs[ed["hom_contig"]][ed["hom_start"] : ed["hom_end"]]
                s = s[:a0] + other + s[b0:]
                ed["delta"] = len(other) - (b0 - a0)
            elif ed["kind"] == "extra":
                extras[f"{ed['id']}_extra"] = seg
                ed["delta"] = 0
        new[c] = s
        acc, table = 0, []
        for ed in sorted(eds, key=lambda x: x["start"]):
            acc += ed["delta"]
            table.append((ed["end"], acc))
        shifts[c] = table

    def to_new(contig, pos):
        """Edited coordinate of an original position that lies outside every edit."""
        acc = 0
        for end, cum in shifts.get(contig, []):
            if pos >= end:
                acc = cum
        return pos + acc

    expected = []
    for ed in edits:
        c, s, e, L = ed["contig"], ed["start"], ed["end"], ed["end"] - ed["start"]
        ns = to_new(c, s - 1) + 1  # edited coordinate of the edit's first base
        hs, he = to_new(ed["hom_contig"], ed["hom_start"]), to_new(ed["hom_contig"], ed["hom_end"] - 1) + 1
        ed["new_start"] = ns
        tag = f"{ed['id']}|{ed['kind']}|{ed['size']}|{ed['stratum']}"
        if ed["kind"] == "dup":
            ed["new_end"] = ns + 2 * L
            expected += [(c, ns, ns + 2 * L, "EXPANSION", tag), (ed["hom_contig"], hs, he, "EXPANSION_PARTNER", tag)]
        elif ed["kind"] == "extra":
            ed["new_end"] = ns + L
            expected += [(f"{ed['id']}_extra", 0, L, "EXPANSION", tag), (c, ns, ns + L, "EXPANSION", tag),
                         (ed["hom_contig"], hs, he, "EXPANSION_PARTNER", tag)]
        elif ed["kind"] == "del":
            ed["new_end"] = ns
            expected += [(ed["hom_contig"], hs, he, "COLLAPSE", tag), (ed["hom_contig"], hs, he, "MISSING_HAPLOTYPE", tag),
                         (c, ns - 1, ns + 1, "FALSE_JOIN", tag)]
        elif ed["kind"] == "swap":
            ed["new_end"] = ns + (ed["hom_end"] - ed["hom_start"])
            expected += [(c, ns, ed["new_end"], "MISSING_HAPLOTYPE", tag), (c, ns, ed["new_end"], "PHASE_SWITCH", tag)]
        elif ed["kind"] == "inv":
            ed["new_end"] = ns + L
            expected += [(c, ns - 1, ns + 1, "MISSING_JUNCTION", tag), (c, ns + L - 1, ns + L + 1, "MISSING_JUNCTION", tag)]
        elif ed["kind"] == "invgap":
            ed["new_end"] = ns + L + 200
            expected += [(c, ns, ns + 100, "MISSING_JUNCTION", tag), (c, ns + 100 + L, ns + 200 + L, "MISSING_JUNCTION", tag)]

    with open(f"{a.out}.fa", "wb") as f:
        for name, s in new.items():
            f.write(b">" + name.encode() + b"\n" + s + b"\n")
        for name, s in extras.items():
            f.write(b">" + name.encode() + b"\n" + s + b"\n")
    cols = ["id", "kind", "size", "stratum", "contig", "start", "end", "new_start", "new_end", "hom_contig", "hom_start", "hom_end"]
    with open(f"{a.out}.edits.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for ed in edits:
            f.write("\t".join(str(ed[k]) for k in cols) + "\n")
    with open(f"{a.out}.expected.bed", "w") as f:
        for row in expected:
            f.write("\t".join(map(str, row)) + "\n")
    with open(f"{a.out}.groups.tsv", "w") as f:
        for name in list(new) + list(extras):
            g = "mat" if "MATERNAL" in name else ("pat" if "PATERNAL" in name else "other")
            f.write(f"{name}\t{g}\n")
    print(f"wrote {a.out}.fa with {len(edits)} edits", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
