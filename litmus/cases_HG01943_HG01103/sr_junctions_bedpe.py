import pysam, sys, re, json
"""Short-read junction evidence for every junction in a bedpe (svreview columns: side1/side2 at 12/13) in one CRAM/BAM.
usage: sr_junctions_bedpe.py cram ref.fa junctions.bedpe label [win=1000] [slack=10] [mapq=20]"""
fn, ref, bedpe, label = sys.argv[1:5]
win = int(sys.argv[5]) if len(sys.argv) > 5 else 1000; slack = int(sys.argv[6]) if len(sys.argv) > 6 else 10; mq = int(sys.argv[7]) if len(sys.argv) > 7 else 20
B = pysam.AlignmentFile(fn, reference_filename=ref)
def touches(side, rs0, re_, p): return abs(re_ - p) <= slack if side == 'L' else abs(rs0 + 1 - p) <= slack
def in_piece(side, r, p): return ((not r.is_reverse) and p - win <= r.reference_start and r.reference_end <= p + slack) if side == 'L' else (r.is_reverse and p - slack <= r.reference_start and r.reference_end <= p + win)
for l in open(bedpe):
    if l.startswith('#'): continue
    t = l.rstrip('\n').split('\t'); c1, p1, c2, p2, name, s1, s2 = t[0], int(t[2]), t[3], int(t[5]), t[6], t[11], t[12]
    pairs, splits = set(), set()
    for (c, p, s, oc, op, os_) in [(c1, p1, s1, c2, p2, s2), (c2, p2, s2, c1, p1, s1)]:
        for r in B.fetch(c, max(0, p - win - 200), p + win + 200):
            if r.is_secondary or r.is_duplicate or r.is_qcfail or r.is_unmapped: continue
            if r.is_paired and not r.mate_is_unmapped and r.mapping_quality >= mq and in_piece(s, r, p):
                mpos = r.next_reference_start; mrev = r.mate_is_reverse; mmq = r.get_tag('MQ') if r.has_tag('MQ') else 60
                ok = (os_ == 'L' and not mrev and op - win <= mpos <= op) or (os_ == 'R' and mrev and op - slack <= mpos <= op + win)
                if r.next_reference_name == oc and ok and mmq >= mq: pairs.add(r.query_name)
            if r.has_tag('SA') and touches(s, r.reference_start, r.reference_end, p):
                for x in r.get_tag('SA').rstrip(';').split(';'):
                    rn, pos, strand, cig, smq, nm = x.split(',')
                    L = sum(int(n) for n, o in re.findall(r'(\d+)([MDN=X])', cig)); st0 = int(pos) - 1
                    if rn == oc and touches(os_, st0, st0 + L, op): splits.add(r.query_name)
    print(json.dumps(dict(sample=label, junction=name, pairs=len(pairs), split_reads=len(splits))), flush=True)
