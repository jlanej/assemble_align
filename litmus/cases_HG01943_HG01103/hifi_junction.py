import pysam, sys, collections, json, re
"""Count reads carrying a junction (split alignment: one part ending at end1 (side L) / starting at (side R)) and
reads straight through each end. Usage: hifi_junction.py bam chrom1 pos1 side1 chrom2 pos2 side2 [slack] [mapq] [flank]"""
bam, c1, p1, s1, c2, p2, s2 = sys.argv[1:8]
p1, p2 = int(p1), int(p2)
slack = int(sys.argv[8]) if len(sys.argv) > 8 else 50
mq = int(sys.argv[9]) if len(sys.argv) > 9 else 20
flank = int(sys.argv[10]) if len(sys.argv) > 10 else 1000
B = pysam.AlignmentFile(bam)
def ends_at(chrom, ref_start0, ref_end, pos, side):
    # side L: piece lies left of pos -> alignment ends at pos ; side R: alignment starts at pos
    if side == 'L': return abs(ref_end - pos) <= slack
    return abs((ref_start0 + 1) - pos) <= slack
def sa_list(r):
    if not r.has_tag('SA'): return []
    out = []
    for x in r.get_tag('SA').rstrip(';').split(';'):
        rn, pos, strand, cig, mapq, nm = x.split(',')
        # reference length from cigar
        L = sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', cig))
        out.append((rn, int(pos) - 1, int(pos) - 1 + L, int(mapq)))
    return out
junction, through = {}, {1: set(), 2: set()}
cells = collections.Counter()
for idx, (c, p, s) in enumerate([(c1, p1, s1), (c2, p2, s2)], start=1):
    for r in B.fetch(c, max(0, p - flank - 100), p + flank + 100):
        if r.is_secondary or r.is_unmapped or r.mapping_quality < mq: continue
        # straight through: aligned block covers p-flank .. p+flank and not clipped near p
        if r.reference_start <= p - flank and r.reference_end >= p + flank:
            through[idx].add(r.query_name)
        if ends_at(c, r.reference_start, r.reference_end, p, s):
            other = (c2, p2, s2) if idx == 1 else (c1, p1, s1)
            for (rn, a, b, m) in sa_list(r):
                if rn == other[0] and m >= mq and ends_at(rn, a, b, other[1], other[2]):
                    junction[r.query_name] = r.get_tag('RG') if r.has_tag('RG') else '.'
res = dict(junction_reads=len(junction), through_1=len(through[1]), through_2=len(through[2]),
           junction_fraction=round(len(junction) / max(1, len(junction) + len(through[1]) + len(through[2])), 4),
           read_groups=dict(collections.Counter(junction.values())))
print(json.dumps(res))
with open(sys.argv[11] if len(sys.argv) > 11 else '/dev/null', 'w') as f:
    for q in sorted(junction): f.write(q + '\n')
