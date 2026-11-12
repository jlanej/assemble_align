#!/usr/bin/env python3
"""Nested random subsampling of a read set to target coverages (each lower coverage is a subset of the
higher ones).  Every read gets a key = blake2b(seed + read name); reads are taken in key order until the
target number of bases (coverage x genome size) is reached.  Two passes, O(#reads) memory.
Usage: 05_nested_subsample.py <reads.fastq.gz> <genome_size_bp> <out_prefix> <cov1,cov2,...|full> [--seed S] [--min-len L]
Writes <out_prefix>.<cov>x.fastq.gz and <out_prefix>.subsample.tsv (target, achieved bases, reads)."""
import sys, gzip, hashlib, argparse
ap = argparse.ArgumentParser(); ap.add_argument('fq'); ap.add_argument('gsize', type=float); ap.add_argument('out')
ap.add_argument('covs'); ap.add_argument('--seed', default='L4-2026'); ap.add_argument('--min-len', type=int, default=0)
a = ap.parse_args()
def key(name): return int.from_bytes(hashlib.blake2b((a.seed + name).encode(), digest_size=8).digest(), 'big')
def records(fn):
    with gzip.open(fn, 'rt') as f:
        while True:
            h = f.readline()
            if not h: return
            s = f.readline(); p = f.readline(); q = f.readline()
            yield h, s, p, q
keys = []
for h, s, p, q in records(a.fq):
    L = len(s) - 1
    if L >= a.min_len: keys.append((key(h[1:].split()[0]), L))
keys.sort(); total = sum(L for _, L in keys)
covs = [c for c in a.covs.split(',')]
thr = {}
for c in covs:
    target = total if c == 'full' else float(c) * a.gsize
    acc = 0; t = None; n = 0
    for k, L in keys:
        if acc >= target: break
        acc += L; t = k; n += 1
    if acc < target and c != 'full':
        print(f'# coverage {c}x not reachable: only {total/a.gsize:.2f}x available', file=sys.stderr); continue
    thr[c] = (t, acc, n)
outs = {c: gzip.open(f'{a.out}.{c}x.fastq.gz', 'wt', compresslevel=1) for c in thr}
for h, s, p, q in records(a.fq):
    if len(s) - 1 < a.min_len: continue
    k = key(h[1:].split()[0])
    for c, (t, _, _) in thr.items():
        if k <= t: outs[c].write(h + s + p + q)
for o in outs.values(): o.close()
import os
new = not os.path.exists(f'{a.out}.subsample.tsv')
with open(f'{a.out}.subsample.tsv', 'a') as f:
    if new: f.write('target\tbases\treads\tachieved_cov\tavailable_cov\n')
    for c, (t, acc, n) in thr.items():
        f.write(f'{c}\t{acc}\t{n}\t{acc/a.gsize:.2f}\t{total/a.gsize:.2f}\n')
print(open(f'{a.out}.subsample.tsv').read())
