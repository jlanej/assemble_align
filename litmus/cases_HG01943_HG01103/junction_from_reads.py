import pysam, sys, re, json, collections
"""For each junction (bedpe: c1 s1 e1 c2 s2 e2 name ... side1 side2 at cols 12,13 of the svreview file), count HiFi reads
split across it, reads straight through each end, and build a read-derived junction sequence (majority string of the
read bases spanning the break, end1 part first). usage: junction_from_reads.py hifi.bam junctions.bedpe out.tsv out_kmers.tsv [k=31]"""
bam, bedpe, out_tsv, out_k = sys.argv[1:5]; K = int(sys.argv[5]) if len(sys.argv) > 5 else 31
B = pysam.AlignmentFile(bam)
comp = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(comp)[::-1]
def cig_parse(c):
    ops = re.findall(r'(\d+)([MIDNSHP=X])', c); ops = [(int(n), o) for n, o in ops]
    lead = ops[0][0] if ops[0][1] in 'SH' else 0; trail = ops[-1][0] if ops[-1][1] in 'SH' else 0
    alen = sum(n for n, o in ops if o in 'MI=X'); rlen = sum(n for n, o in ops if o in 'MDN=X')
    return lead, alen, trail, rlen
def touches(side, rs1, re_, p, slack=50):
    return abs(re_ - p) <= slack if side == 'L' else abs(rs1 - p) <= slack
J = []
for l in open(bedpe):
    if l.startswith('#'): continue
    t = l.rstrip('\n').split('\t')
    J.append(dict(c1=t[0], p1=int(t[2]), c2=t[3], p2=int(t[5]), name=t[6], s1=t[11], s2=t[12], repo_hifi=int(t[14])))
rows = []; kmers = []
for j in J:
    junc = {}; thr = {1: set(), 2: set()}
    for idx, (c, p, s) in enumerate([(j['c1'], j['p1'], j['s1']), (j['c2'], j['p2'], j['s2'])], 1):
        for r in B.fetch(c, max(0, p - 1200), p + 1200):
            if r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
            if r.reference_start <= p - 1000 and r.reference_end >= p + 1000: thr[idx].add(r.query_name)
    # junction reads: primary records with an SA part such that one part touches end1 and the other touches end2
    seen = set()
    for c, p in [(j['c1'], j['p1']), (j['c2'], j['p2'])]:
        for r in B.fetch(c, max(0, p - 30000), p + 30000):
            if r.is_secondary or r.is_supplementary or r.query_name in seen or not r.has_tag('SA') or r.mapping_quality < 20: continue
            full = r.query_sequence
            if full is None: continue
            qlen = len(full)
            parts = []
            lead, alen, trail, rlen = cig_parse(r.cigarstring)
            parts.append((r.reference_name, r.reference_start + 1, r.reference_end, r.is_reverse, lead, alen, r.mapping_quality))
            for x in r.get_tag('SA').rstrip(';').split(';'):
                rn, pos, strand, cg, mq, nm = x.split(',')
                lead2, alen2, trail2, rlen2 = cig_parse(cg)
                parts.append((rn, int(pos), int(pos) + rlen2 - 1, strand == '-', lead2, alen2, int(mq)))
            def qint(pt):
                rn, rs, re_, rev, lead_, alen_, mq = pt
                return (qlen - lead_ - alen_, qlen - lead_) if rev else (lead_, lead_ + alen_)
            a = [pt for pt in parts if pt[0] == j['c1'] and pt[6] >= 20 and touches(j['s1'], pt[1], pt[2], j['p1'])]
            b = [pt for pt in parts if pt[0] == j['c2'] and pt[6] >= 20 and touches(j['s2'], pt[1], pt[2], j['p2'])]
            hit = None
            for pa in a:
                for pb in b:
                    if pa is pb: continue
                    qa, qb = qint(pa), qint(pb)
                    if qa[1] <= qb[0] + 200 and qa[1] >= qb[0] - 200: hit = (qa, qb, False)
                    elif qb[1] <= qa[0] + 200 and qb[1] >= qa[0] - 200: hit = (qa, qb, True)
            if hit:
                seen.add(r.query_name)
                orig = rc(full) if r.is_reverse else full
                qa, qb, swapped = hit
                if swapped:   # end2 part first in read: flip read so end1 part comes first
                    orig = rc(orig); qa = (qlen - qa[1], qlen - qa[0]); qb = (qlen - qb[1], qlen - qb[0])
                lo, hi = min(qa[1], qb[0]), max(qa[1], qb[0])
                junc[r.query_name] = (orig[max(0, lo - 30): hi + 30], hi - lo, qa[1] - qb[0])
    T1, T2, Jn = len(thr[1] - set(junc)), len(thr[2] - set(junc)), len(junc)
    f = 2 * Jn / (Jn + (T1 + T2) / 2) if Jn + T1 + T2 else float('nan')
    # consensus: most common break-window string
    cnt = collections.Counter(v[0] for v in junc.values() if len(v[0]) >= 60)
    cons, ncons = (cnt.most_common(1)[0] if cnt else ('', 0))
    gap = collections.Counter(v[2] for v in junc.values()).most_common(1)[0][0] if junc else 0
    rows.append(dict(junction=j['name'], hifi_junction_reads=Jn, repo_hifi_reads=j['repo_hifi'], through_end1=T1, through_end2=T2,
                     f_hifi=round(f, 3), consensus_support=ncons, overlap_or_gap_bp=gap, consensus=cons))
    if cons:
        L = len(cons); mid0 = 30; mid1 = L - 30   # break region [mid0, mid1)
        for st in range(max(0, mid1 + 3 - K), min(mid0 - 3, L - K) + 1):
            kmers.append((f"{j['name']}_{st}", cons[st:st + K]))
with open(out_tsv, 'w') as f:
    ks = list(rows[0].keys()); f.write('\t'.join(ks) + '\n')
    for r in rows: f.write('\t'.join(str(r[k]) for k in ks) + '\n')
with open(out_k, 'w') as f:
    for n, s in kmers: f.write(f'{n}\t{s}\n')
for r in rows: print({k: v for k, v in r.items() if k != 'consensus'})
print('kmers', len(kmers))
