import pysam, sys, re, json
"""Long reads (CRAM/BAM) at every junction; a junction with 0 fetched records is reported as NO_READS_FETCHED (EBI range reads can truncate silently) of a bedpe: reads split across the junction (SA, both parts MAPQ>=20, ends
within 100 bp) and reads straight through each end (+-1 kb). usage: lr_junctions_bedpe.py aln ref.fa junctions.bedpe label"""
fn, ref, bedpe, label = sys.argv[1:5]
B = pysam.AlignmentFile(fn, reference_filename=ref)
def rlen(c): return sum(int(n) for n, o in re.findall(r'(\d+)([MDN=X])', c))
def touches(side, rs0, re_, p, slack=100): return abs(re_ - p) <= slack if side == 'L' else abs(rs0 + 1 - p) <= slack
for l in open(bedpe):
    if l.startswith('#'): continue
    t = l.rstrip('\n').split('\t'); c1, p1, c2, p2, name, s1, s2 = t[0], int(t[2]), t[3], int(t[5]), t[6], t[11], t[12]
    junc, thr = set(), {1: set(), 2: set()}; fetched = 0
    for e, (c, p, s, oc, op, os_) in enumerate([(c1, p1, s1, c2, p2, s2), (c2, p2, s2, c1, p1, s1)], 1):
        for r in B.fetch(c, max(0, p - 1100), p + 1100):
            fetched += 1
            if r.is_secondary or r.is_unmapped or r.mapping_quality < 20: continue
            if r.reference_start <= p - 1000 and r.reference_end >= p + 1000 and not r.is_supplementary: thr[e].add(r.query_name)
            if r.has_tag('SA') and touches(s, r.reference_start, r.reference_end, p):
                for x in r.get_tag('SA').rstrip(';').split(';'):
                    rn, pos, strand, cg, smq, nm = x.split(',')
                    if rn == oc and int(smq) >= 20 and touches(os_, int(pos) - 1, int(pos) - 1 + rlen(cg), op): junc.add(r.query_name)
    print(json.dumps(dict(sample=label, junction=name, junction_reads=len(junc), through_end1=len(thr[1]), through_end2=len(thr[2]), records_fetched=fetched,
                     status=('NO_READS_FETCHED (truncated remote read? retry)' if fetched == 0 else 'ok'))), flush=True)
