#!/usr/bin/env python3
"""Copy number from CRAM indexes alone (no reads, no reference): each NYGC CRAM slice holds 10,000 reads of 150 bp,
so mean depth over a region = 1.5e6 * n_slices / sum(slice spans). The child's segment depth is normalised by its own
flank (or a chosen normaliser region) and by the parents' same ratio. usage:
  crai_cn.py --child URL.crai --father URL.crai --mother URL.crai --regions regions.tsv [--norm chrom:start-end] --out out.tsv
regions.tsv: name<TAB>chrom<TAB>start<TAB>end (1-based; 0-based BED is accepted with --bed). chrom index order of the
1000G GRCh38 analysis set: chr1..chr22 = 0..21, chrX = 22."""
import argparse, gzip, io, urllib.request
CH = {f'chr{i}': i - 1 for i in range(1, 23)}; CH['chrX'] = 22
def load(u):
    raw = urllib.request.urlopen(u, timeout=300).read() if u.startswith('http') else open(u, 'rb').read()
    rows = []
    for line in gzip.decompress(raw).decode().splitlines():
        t = line.split('\t'); rows.append((int(t[0]), int(t[1]), int(t[2])))
    return rows
def depth(rows, chrom, lo, hi):
    ref = CH[chrom]; sel = [(s, sp) for r, s, sp in rows if r == ref and sp > 0 and s >= lo and s + sp <= hi]
    return (1.5e6 * len(sel) / sum(sp for _, sp in sel), len(sel)) if sel else (float('nan'), 0)
ap = argparse.ArgumentParser(); ap.add_argument('--child'); ap.add_argument('--father'); ap.add_argument('--mother')
ap.add_argument('--regions', required=True); ap.add_argument('--bed', action='store_true'); ap.add_argument('--norm', default='chr2:20000000-140000000'); ap.add_argument('--out', required=True)
a = ap.parse_args()
C = {k: load(getattr(a, k)) for k in ('child', 'father', 'mother')}
nc, ns = a.norm.split(':'); n0, n1 = map(int, ns.split('-'))
norm = {k: depth(v, nc, n0, n1)[0] for k, v in C.items()}
with open(a.out, 'w') as f:
    f.write('name\tchrom\tstart\tend\tchild_depth\tfather_depth\tmother_depth\tchild_cn\tfather_cn\tmother_cn\tn_slices_child\n')
    for line in open(a.regions):
        if line.startswith('#') or not line.strip(): continue
        t = line.rstrip('\n').split('\t')
        if a.bed: chrom, s, e, name = t[0], int(t[1]) + 1, int(t[2]), (t[3] if len(t) > 3 else f'{t[0]}:{t[1]}-{t[2]}')
        else: name, chrom, s, e = t[0], t[1], int(t[2]), int(t[3])
        if e - s < 100000: continue
        d = {k: depth(v, chrom, s + 20000, e - 20000) for k, v in C.items()}
        r = {k: d[k][0] / norm[k] for k in d}
        cn = 2 * r['child'] / ((r['father'] + r['mother']) / 2)
        f.write(f"{name}\t{chrom}\t{s}\t{e}\t{d['child'][0]:.2f}\t{d['father'][0]:.2f}\t{d['mother'][0]:.2f}\t{cn:.3f}\t{2*r['father']:.3f}\t{2*r['mother']:.3f}\t{d['child'][1]}\n")
