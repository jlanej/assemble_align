import pysam, sys, re, collections, json
"""Read-to-assembly QC: per-bin depth (all, MAPQ>=20, MAPQ0 share), depth by read parent, long soft-clip ends, reads
spanning each N gap. usage: asm_read_qc.py reads_on_asm.bam asm.fa read_parent.tsv out_prefix [bin=5000] [clip=500]"""
bam, fa, rp_tsv, outp = sys.argv[1:5]
bin_ = int(sys.argv[5]) if len(sys.argv) > 5 else 5000
clipmin = int(sys.argv[6]) if len(sys.argv) > 6 else 500
par = {l.split()[0]: l.split()[3] for l in open(rp_tsv)}
F = pysam.FastaFile(fa); B = pysam.AlignmentFile(bam)
rows = []; gaps = []
for ctg in F.references:
    L = F.get_reference_length(ctg); seq = F.fetch(ctg).upper()
    nb = (L + bin_ - 1) // bin_
    d_all = [0]*nb; d_q20 = [0]*nb; d_q0 = [0]*nb; d_mat = [0]*nb; d_pat = [0]*nb; clips = [0]*nb
    for r in B.fetch(ctg):
        if r.is_secondary or r.is_unmapped: continue
        cls = par.get(r.query_name, 'unph')
        for (s, e) in r.get_blocks():
            for b in range(s // bin_, (e - 1) // bin_ + 1):
                ov = min(e, (b + 1) * bin_) - max(s, b * bin_)
                d_all[b] += ov
                if r.mapping_quality >= 20: d_q20[b] += ov
                if r.mapping_quality == 0: d_q0[b] += ov
                if cls == 'mat': d_mat[b] += ov
                if cls == 'pat': d_pat[b] += ov
        ct = r.cigartuples
        if ct[0][0] in (4, 5) and ct[0][1] >= clipmin: clips[r.reference_start // bin_] += 1
        if ct[-1][0] in (4, 5) and ct[-1][1] >= clipmin: clips[(r.reference_end - 1) // bin_] += 1
    for b in range(nb):
        w = min(L, (b + 1) * bin_) - b * bin_
        rows.append((ctg, b * bin_, b * bin_ + w, d_all[b] / w, d_q20[b] / w, d_q0[b] / max(1, d_all[b]), d_mat[b] / w, d_pat[b] / w, clips[b]))
    for m in re.finditer('N+', seq):
        gs, ge = m.start(), m.end()
        span = sum(1 for r in B.fetch(ctg, max(0, gs - 1000), ge + 1000) if not r.is_secondary and r.reference_start <= gs - 1000 and r.reference_end >= ge + 1000)
        # base-level depth in 2 kb either side of the gap
        def mean_depth(a, b):
            tot = 0
            for col in B.pileup(ctg, a, b, truncate=True, stepper='nofilter', min_mapping_quality=0): tot += col.nsegments
            return tot / (b - a)
        gaps.append(dict(contig=ctg, gap_start=gs + 1, gap_end=ge, reads_spanning_gap_1kb=span,
                         depth_left_2kb=round(mean_depth(max(0, gs - 2000), gs), 2), depth_right_2kb=round(mean_depth(ge, ge + 2000), 2)))
with open(outp + '.bins.tsv', 'w') as f:
    f.write('contig\tstart\tend\tdepth\tdepth_mq20\tfrac_mq0\tdepth_mat\tdepth_pat\tclip_ends\n')
    for r in rows: f.write('\t'.join(str(round(x, 3)) if isinstance(x, float) else str(x) for x in r) + '\n')
json.dump(gaps, open(outp + '.gaps.json', 'w'), indent=1)
print(json.dumps(gaps, indent=1))
