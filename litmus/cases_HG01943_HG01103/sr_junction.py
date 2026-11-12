import pysam, sys, json, re
"""Short-read evidence for a junction joining end1 (chrom,pos,side) to end2. side L: piece lies left of pos (reads
point right, + strand); side R: piece lies right of pos (reads point left, - strand).
usage: sr_junction.py BAM/CRAM c1 p1 s1 c2 p2 s2 [ref.fa] [window=1000] [slack=10] [mapq=20]"""
a = sys.argv
fn, c1, p1, s1, c2, p2, s2 = a[1], a[2], int(a[3]), a[4], a[5], int(a[6]), a[7]
ref = a[8] if len(a) > 8 and a[8] != '-' else None
win = int(a[9]) if len(a) > 9 else 1000
slack = int(a[10]) if len(a) > 10 else 10
mq = int(a[11]) if len(a) > 11 else 20
B = pysam.AlignmentFile(fn, reference_filename=ref)
def in_piece(chrom, pos, side, rs, re_, rev, c, p):
    if chrom != c: return False
    if side == 'L': return (not rev) and p - win <= rs and re_ <= p + slack
    return rev and p - slack <= rs and re_ <= p + win
def touches(side, rs0, re_, p):
    return abs(re_ - p) <= slack if side == 'L' else abs(rs0 + 1 - p) <= slack
pairs, splits, splits_any = set(), set(), set()
for (c, p, s, oc, op, os_) in [(c1, p1, s1, c2, p2, s2), (c2, p2, s2, c1, p1, s1)]:
    for r in B.fetch(c, max(0, p - win - 200), p + win + 200):
        if r.is_secondary or r.is_duplicate or r.is_qcfail or r.is_unmapped: continue
        # discordant pair
        if r.is_paired and not r.mate_is_unmapped and r.mapping_quality >= mq:
            if in_piece(c, p, s, r.reference_start, r.reference_end, r.is_reverse, c, p):
                mpos = r.next_reference_start
                mrev = r.mate_is_reverse
                mchrom = r.next_reference_name
                mmq = r.get_tag('MQ') if r.has_tag('MQ') else 60
                ok_side = (os_ == 'L' and (not mrev) and op - win <= mpos <= op) or (os_ == 'R' and mrev and op - slack <= mpos <= op + win)
                if mchrom == oc and ok_side and mmq >= mq:
                    pairs.add(r.query_name)
        # split read
        if r.has_tag('SA') and touches(s, r.reference_start, r.reference_end, p):
            for x in r.get_tag('SA').rstrip(';').split(';'):
                rn, pos, strand, cig, smq, nm = x.split(',')
                L = sum(int(n) for n, o in re.findall(r'(\d+)([MDN=X])', cig))
                st0 = int(pos) - 1
                if rn == oc and touches(os_, st0, st0 + L, op):
                    splits_any.add(r.query_name)
                    if int(smq) >= mq and r.mapping_quality >= mq: splits.add(r.query_name)
print(json.dumps(dict(file=fn.split('/')[-1], pairs=len(pairs), split_reads_mq20=len(splits), split_reads_anymq=len(splits_any))))
