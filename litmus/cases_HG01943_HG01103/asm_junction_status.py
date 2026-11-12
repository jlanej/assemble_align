import pysam, sys, re, json, collections
"""Classify assembly-to-reference alignments at each junction end: junction held (same contig, adjacent query intervals,
one part ending at each end), straight through (record spans end +-1 kb), or contig/alignment end.
usage: asm_junction_status.py junctions.bedpe segments.bed out.tsv hap1.bam [hap2.bam ...]"""
bedpe, segbed, out = sys.argv[1:4]; bams = sys.argv[4:]
def qint(r):
    ct = r.cigartuples; lead = ct[0][1] if ct[0][0] in (4, 5) else 0; trail = ct[-1][1] if ct[-1][0] in (4, 5) else 0
    alen = sum(n for o, n in ct if o in (0, 1, 7, 8)); qlen = lead + alen + trail
    return (qlen - lead - alen, qlen - lead) if r.is_reverse else (lead, lead + alen)
def touches(side, r, p, slack=100):
    return abs(r.reference_end - p) <= slack if side == 'L' else abs(r.reference_start + 1 - p) <= slack
J = []
for l in open(bedpe):
    if l.startswith('#'): continue
    t = l.rstrip('\n').split('\t'); J.append(dict(c1=t[0], p1=int(t[2]), c2=t[3], p2=int(t[5]), name=t[6], s1=t[11], s2=t[12], repo=(t[19] if len(t) > 19 else '-')))
rows = []
for fn in bams:
    hap = fn.split('/')[-1].replace('.bam', '')
    B = pysam.AlignmentFile(fn)
    for j in J:
        recs = {}
        for e, (c, p, s) in enumerate([(j['c1'], j['p1'], j['s1']), (j['c2'], j['p2'], j['s2'])], 1):
            recs[e] = [r for r in B.fetch(c, max(0, p - 1000), p + 1000) if not r.is_secondary and not r.is_unmapped]
        thr = {e: sorted({r.query_name for r in recs[e] if r.reference_start <= p - 1000 and r.reference_end >= p + 1000})
               for e, p in [(1, j['p1']), (2, j['p2'])]}
        ends = {e: [r for r in recs[e] if touches(s, r, p)] for e, (p, s) in [(1, (j['p1'], j['s1'])), (2, (j['p2'], j['s2']))]}
        held = set()
        for a in ends[1]:
            for b in ends[2]:
                if a.query_name != b.query_name or (a.reference_start == b.reference_start and a.reference_name == b.reference_name): continue
                qa, qb = qint(a), qint(b)
                if min(abs(qa[1] - qb[0]), abs(qb[1] - qa[0])) <= 200: held.add(a.query_name)
        status = 'junction_held' if held else ('through_both' if thr[1] and thr[2] else 'through_one' if (thr[1] or thr[2]) else 'contig_ends')
        rows.append(dict(haplotype=hap, junction=j['name'], status=status, contig_holding=','.join(sorted(held)) or '-',
                         through_end1=','.join(thr[1]) or '-', through_end2=','.join(thr[2]) or '-',
                         ending_end1=','.join(sorted({r.query_name for r in ends[1]})) or '-', ending_end2=','.join(sorted({r.query_name for r in ends[2]})) or '-',
                         repo_assembly_hap1=j['repo']))
# segment copy numbers
seg_rows = []
for l in open(segbed):
    if l.startswith('#'): continue
    c, s, e, piece = l.split('\t')[:4]; s, e = int(s), int(e)
    if e - s < 100000: continue
    cn = {}
    for fn in bams:
        hap = fn.split('/')[-1].replace('.bam', ''); B = pysam.AlignmentFile(fn); cov = 0
        for r in B.fetch(c, s, e):
            if r.is_secondary or r.is_unmapped: continue
            for a, b in r.get_blocks():
                ov = min(b, e) - max(a, s)
                if ov > 0: cov += ov
        cn[hap] = round(cov / (e - s), 3)
    seg_rows.append(dict(chrom=c, start=s, end=e, piece=piece, **cn))
with open(out, 'w') as f:
    ks = list(rows[0].keys()); f.write('\t'.join(ks) + '\n')
    for r in rows: f.write('\t'.join(str(r[k]) for k in ks) + '\n')
with open(out.replace('.tsv', '.segments.tsv'), 'w') as f:
    ks = list(seg_rows[0].keys()); f.write('\t'.join(ks) + '\n')
    for r in seg_rows: f.write('\t'.join(str(r[k]) for k in ks) + '\n')
for r in rows: print(r['haplotype'], r['junction'], r['status'], r['contig_holding'], 'thr1', r['through_end1'], 'thr2', r['through_end2'], 'repo', r['repo_assembly_hap1'])
for r in seg_rows: print(r)
