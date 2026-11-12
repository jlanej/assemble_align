#!/usr/bin/env python3
"""Per-bin k-mer copy number from meryl-lookup -wig-count tracks (one wig per k-mer database, all on the same sequence):
CN_reads = median read count of k-mers present in both reads and assembly / haploid k-mer depth; CN_asm = median count in
the assembly database. usage: kmer_cn_track.py --offset 46400000 --haploid-depth 15 --bin 10000 --out bins.tsv child.wig father.wig mother.wig asm.wig"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--offset', type=int, default=0); ap.add_argument('--haploid-depth', type=float, required=True)
ap.add_argument('--bin', type=int, default=10000); ap.add_argument('--min-kmers', type=int, default=2000); ap.add_argument('--out', required=True); ap.add_argument('wigs', nargs=4)
a = ap.parse_args()
def load(f):
    vals = {}
    for line in open(f):
        if line[0] in 'tvf': continue
        p, v = line.split(); vals[int(p)] = int(float(v))
    L = max(vals) + 1; arr = np.zeros(L, dtype=np.int32)
    for p, v in vals.items(): arr[p - 1] = v
    return arr
c, fa, mo, asm = (load(w) for w in a.wigs); L = min(map(len, (c, fa, mo, asm)))
with open(a.out, 'w') as f:
    f.write('start\tend\tn_kmers\tcn_child_reads\tcn_father_reads\tcn_mother_reads\tcn_asm\n')
    for b in range(0, L, a.bin):
        sl = slice(b, min(L, b + a.bin)); m = (asm[sl] > 0) & (c[sl] > 0)
        if m.sum() < a.min_kmers: continue
        h = a.haploid_depth
        f.write(f"{a.offset+b+1}\t{a.offset+b+a.bin}\t{int(m.sum())}\t{np.median(c[sl][m])/h:.2f}\t{np.median(fa[sl][m])/h:.2f}\t{np.median(mo[sl][m])/h:.2f}\t{np.median(asm[sl][m]):.1f}\n")
