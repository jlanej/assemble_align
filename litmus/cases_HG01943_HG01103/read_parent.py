import pysam, sys
"""Classify every primary read of a GRCh38 BAM by trio alleles. sites: POS REF ALT GTchild GTfather GTmother"""
sites, bam, chrom, out = sys.argv[1:5]
monly, fonly = {}, {}
for l in open(sites):
    pos, ref, alt, gc, gf, gm = l.split()
    if ',' in alt: continue
    fa = set(gf.replace('|','/').split('/')); ma = set(gm.replace('|','/').split('/'))
    if '.' in fa or '.' in ma: continue
    al = {'0': ref, '1': alt}
    for a_ in ('0','1'):
        if a_ in ma and a_ not in fa: monly[int(pos)] = al[a_]
        if a_ in fa and a_ not in ma: fonly[int(pos)] = al[a_]
B = pysam.AlignmentFile(bam); res = {}
for r in B.fetch(chrom):
    if r.is_secondary or r.is_unmapped: continue
    m, f = res.get(r.query_name, (0, 0))
    q = r.query_sequence
    if q is None: continue
    for qp, rp in r.get_aligned_pairs(matches_only=True):
        p = rp + 1
        if p in monly and q[qp] == monly[p]: m += 1
        if p in fonly and q[qp] == fonly[p]: f += 1
    res[r.query_name] = (m, f)
with open(out, 'w') as fo:
    for k, (m, f) in res.items():
        cls = 'mat' if m >= 2 and f <= 0.05 * m else 'pat' if f >= 2 and m <= 0.05 * f else 'unph'
        fo.write(f'{k}\t{m}\t{f}\t{cls}\n')
