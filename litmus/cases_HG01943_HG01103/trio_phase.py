import pysam, sys, json, collections
"""Parent-of-origin from trio SNVs. Sites: bcftools query POS REF ALT GTchild GTfather GTmother (QUAL>=30).
Outputs: child maternal allele share at parent-opposite sites (in/out of a segment) from short-read AD,
parent assignment of listed reads (HiFi BAM), and parent alleles carried by each assembly piece (contig BAM)."""
sites_tsv, sr_child, hifi_bam, junc_list, asm_bams, chrom, seg_s, seg_e = sys.argv[1:9]
seg_s, seg_e = int(seg_s), int(seg_e)
po, fonly, monly = {}, {}, {}
for l in open(sites_tsv):
    pos, ref, alt, gc, gf, gm = l.split()
    if ',' in alt: continue
    pos = int(pos)
    fa = set(gf.replace('|','/').split('/')); ma = set(gm.replace('|','/').split('/'))
    if '.' in fa or '.' in ma: continue
    al = {'0': ref, '1': alt}
    # parent-opposite: both parents homozygous for different alleles
    if len(fa) == 1 and len(ma) == 1 and fa != ma and gc in ('0/1','0|1','1|0'):
        po[pos] = (al[list(ma)[0]], al[list(fa)[0]])  # (maternal, paternal) base
    # alleles private to one parent
    for a_ in ('0','1'):
        if a_ in ma and a_ not in fa: monly[pos] = al[a_]
        if a_ in fa and a_ not in ma: fonly[pos] = al[a_]
# 1. child maternal allele share (short reads) at PO sites
B = pysam.AlignmentFile(sr_child)
inn, out = [0,0], [0,0]
for col in B.pileup(chrom, min(po)-1, max(po)+1, truncate=True, min_base_quality=20, min_mapping_quality=20, stepper='nofilter'):
    p = col.reference_pos + 1
    if p not in po: continue
    m, f = po[p]; nm = nf = 0
    for pr in col.pileups:
        if pr.is_del or pr.is_refskip or pr.alignment.is_duplicate or pr.alignment.is_secondary: continue
        b = pr.alignment.query_sequence[pr.query_position]
        nm += b == m; nf += b == f
    tgt = inn if seg_s + 5000 <= p <= seg_e - 5000 else (out if (p < seg_s - 5000 or p > seg_e + 5000) else None)
    if tgt is not None: tgt[0] += nm; tgt[1] += nf
res = {'po_sites': len(po), 'mother_only_sites': len(monly), 'father_only_sites': len(fonly),
       'child_maternal_share_inside': round(inn[0]/max(1,sum(inn)),3), 'child_maternal_share_outside': round(out[0]/max(1,sum(out)),3)}
# 2. HiFi junction reads parent
junc = set(x.strip() for x in open(junc_list)) if junc_list != '-' else set()
H = pysam.AlignmentFile(hifi_bam)
def read_parent(r):
    nm = nf = 0
    for qp, rp in r.get_aligned_pairs(matches_only=True):
        p = rp + 1
        if p in monly or p in fonly:
            b = r.query_sequence[qp]
            if p in monly and b == monly[p]: nm += 1
            if p in fonly and b == fonly[p]: nf += 1
    return nm, nf
jp = collections.Counter(); seen = set()
for r in H.fetch(chrom, seg_s - 40000, seg_e + 40000):
    if r.query_name in junc and not r.is_secondary and r.query_name + str(r.reference_start) not in seen:
        seen.add(r.query_name + str(r.reference_start))
        nm, nf = read_parent(r)
        jp[r.query_name] = (jp[r.query_name][0] + nm, jp[r.query_name][1] + nf) if r.query_name in jp else (nm, nf)
cls = collections.Counter('mother' if m >= 2 and f == 0 else 'father' if f >= 2 and m == 0 else 'unphased' for m, f in jp.values())
res['junction_reads_parent'] = dict(cls)
# 3. assembly pieces
pieces = []
for fn in asm_bams.split(','):
    A = pysam.AlignmentFile(fn)
    for r in A.fetch(chrom, seg_s - 1, seg_e):
        nm, nf = 0, 0
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            p = rp + 1
            if not (seg_s <= p <= seg_e): continue
            if p in monly or p in fonly:
                b = r.query_sequence[qp]
                nm += p in monly and b == monly[p]; nf += p in fonly and b == fonly[p]
        pieces.append(dict(contig=r.query_name, ref_start=r.reference_start + 1, ref_end=r.reference_end,
                           mother_only_alleles=nm, father_only_alleles=nf))
res['pieces'] = pieces
print(json.dumps(res, indent=1))
