#!/usr/bin/env python3
"""L2/L3 pilot: residual long-read errors versus the sample's own assembly, by sequence
context, and empirical-vs-reported base-quality calibration.

  l2l3_lrerr.py --pair chosen_h1.bam hap1.fa CONTIG1[,CONTIG..] \
                --pair chosen_h2.bam hap2.fa CONTIG2[,..] --minmapq 20 --out prefix

BAMs: minimap2 --eqx alignments (records of reads assigned to that haplotype), sorted or
unsorted. Only primary+supplementary records with MAPQ>=minmapq on the listed contigs
are used. Context of each assembly base: homopolymer run length (HP) of the run that
contains it; for HP<=3, membership in a 2-6 bp-unit STR (pytrf; >=12 bp and >=6/4/3/3/3
copies for 2/3/4/5/6-mers); otherwise HP1 or HP2-3.
Mismatches and deletions are attributed to the context of the (first) affected assembly
base; an insertion to the adjacent run whose base equals the first inserted base
(left run first), else to the base right of the insertion point.
Outputs <prefix>.context.tsv and <prefix>.qv_calibration.tsv
"""
import argparse, sys
import numpy as np
import pysam
import pytrf

CLASSES = ['HP1', 'HP2-3', 'STR', 'HP4-5', 'HP6-7', 'HP8-10', 'HP11-15', 'HP16+', 'N']


def context_array(seq):
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    n = len(a)
    change = np.empty(n, bool)
    change[0] = True
    change[1:] = a[1:] != a[:-1]
    starts = np.flatnonzero(change)
    lens = np.diff(np.append(starts, n))
    hp = np.repeat(lens, lens)
    strm = np.zeros(n, bool)
    for r in pytrf.STRFinder('x', seq.upper(), 100000, 6, 4, 3, 3, 3):
        if r.length >= 12 and len(r.motif) >= 2:
            strm[r.start - 1:r.end] = True
    cls = np.zeros(n, np.int8)
    cls[(hp >= 2) & (hp <= 3)] = 1
    cls[(hp <= 3) & strm] = 2
    cls[(hp >= 4) & (hp <= 5)] = 3
    cls[(hp >= 6) & (hp <= 7)] = 4
    cls[(hp >= 8) & (hp <= 10)] = 5
    cls[(hp >= 11) & (hp <= 15)] = 6
    cls[hp >= 16] = 7
    cls[a == ord('N')] = 8
    return a, cls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pair', nargs=3, action='append', required=True, metavar=('BAM', 'FASTA', 'CONTIGS'))
    ap.add_argument('--minmapq', type=int, default=20)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    K = len(CLASSES)
    tot = dict(ref=np.zeros(K), mm=np.zeros(K), dele=np.zeros(K), dev=np.zeros(K),
               ins=np.zeros(K), iev=np.zeros(K))
    qv_tot = np.zeros(100, np.int64)
    qv_err = np.zeros(100, np.int64)
    qv_mm = np.zeros(100, np.int64)
    nrec = nbases = 0
    for bam, fasta, contigs in a.pair:
        fa = pysam.FastaFile(fasta)
        ctgs = contigs.split(',')
        ctx = {}
        for c in ctgs:
            ctx[c] = context_array(fa.fetch(c))
        cov = {c: np.zeros(len(ctx[c][0]) + 1, np.int32) for c in ctgs}
        mm_pos = {c: [] for c in ctgs}
        del_ev = {c: [] for c in ctgs}
        ins_ev = {c: [] for c in ctgs}
        b = pysam.AlignmentFile(bam, 'rb', check_sq=False)
        for r in b.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.mapping_quality < a.minmapq:
                continue
            c = r.reference_name
            if c not in ctx:
                continue
            seqa, _ = ctx[c]
            qs = r.query_sequence
            qq = r.query_qualities
            ql = len(qs) if qs else 0
            aligned = np.zeros(ql, bool)
            err = np.zeros(ql, bool)
            mmq = np.zeros(ql, bool)
            rp, qp = r.reference_start, 0
            cv = cov[c]
            for op, n in r.cigartuples:
                if op == 7:      # =
                    cv[rp] += 1
                    cv[rp + n] -= 1
                    aligned[qp:qp + n] = True
                    rp += n
                    qp += n
                elif op == 8:    # X
                    cv[rp] += 1
                    cv[rp + n] -= 1
                    mm_pos[c].extend(range(rp, rp + n))
                    aligned[qp:qp + n] = True
                    err[qp:qp + n] = True
                    mmq[qp:qp + n] = True
                    rp += n
                    qp += n
                elif op == 1:    # I
                    ib = ord(qs[qp])
                    if rp > 0 and seqa[rp - 1] == ib:
                        cpos = rp - 1
                    else:
                        cpos = min(rp, len(seqa) - 1)
                    ins_ev[c].append((cpos, n))
                    aligned[qp:qp + n] = True
                    err[qp:qp + n] = True
                    qp += n
                elif op == 2:    # D
                    cv[rp] += 1
                    cv[rp + n] -= 1
                    del_ev[c].append((rp, n))
                    rp += n
                elif op == 4:
                    qp += n
                elif op == 0:
                    raise SystemExit('CIGAR M found: map with --eqx')
            nrec += 1
            nbases += int(aligned.sum())
            if qq is not None and ql:
                qa = np.asarray(qq)
                qv_tot += np.bincount(np.minimum(qa[aligned], 99), minlength=100)
                qv_err += np.bincount(np.minimum(qa[aligned & err], 99), minlength=100)
                qv_mm += np.bincount(np.minimum(qa[mmq], 99), minlength=100)
        for c in ctgs:
            seqa, cls = ctx[c]
            cover = np.cumsum(cov[c])[:-1]
            tot['ref'] += np.bincount(cls, weights=cover, minlength=K)[:K]
            if mm_pos[c]:
                tot['mm'] += np.bincount(cls[np.array(mm_pos[c])], minlength=K)[:K]
            if del_ev[c]:
                d = np.array(del_ev[c])
                tot['dev'] += np.bincount(cls[d[:, 0]], minlength=K)[:K]
                tot['dele'] += np.bincount(cls[d[:, 0]], weights=d[:, 1], minlength=K)[:K]
            if ins_ev[c]:
                d = np.array(ins_ev[c])
                tot['iev'] += np.bincount(cls[d[:, 0]], minlength=K)[:K]
                tot['ins'] += np.bincount(cls[d[:, 0]], weights=d[:, 1], minlength=K)[:K]
    with open(a.out + '.context.tsv', 'w') as fo:
        fo.write('context\tasm_bases_covered\tmismatch_bases\tins_events\tins_bases\tdel_events\tdel_bases\t'
                 'err_per_kb\tmm_per_kb\tins_ev_per_kb\tdel_ev_per_kb\tfrac_of_bases\n')
        allref = tot['ref'].sum()
        for k, name in enumerate(CLASSES):
            R = tot['ref'][k]
            if R == 0:
                continue
            e = tot['mm'][k] + tot['ins'][k] + tot['dele'][k]
            fo.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\n' % (
                name, R, tot['mm'][k], tot['iev'][k], tot['ins'][k], tot['dev'][k], tot['dele'][k],
                1000 * e / R, 1000 * tot['mm'][k] / R, 1000 * tot['iev'][k] / R, 1000 * tot['dev'][k] / R, R / allref))
        R = allref
        e = tot['mm'].sum() + tot['ins'].sum() + tot['dele'].sum()
        fo.write('ALL\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t1\n' % (
            R, tot['mm'].sum(), tot['iev'].sum(), tot['ins'].sum(), tot['dev'].sum(), tot['dele'].sum(),
            1000 * e / R, 1000 * tot['mm'].sum() / R, 1000 * tot['iev'].sum() / R, 1000 * tot['dev'].sum() / R))
    with open(a.out + '.qv_calibration.tsv', 'w') as fo:
        fo.write('reported_qv\tbases\terrors_mm_or_ins\tmismatches\tempirical_qv\n')
        for q in range(100):
            if qv_tot[q] == 0:
                continue
            e = qv_err[q]
            emp = -10 * np.log10(max(e, 0.5) / qv_tot[q])
            fo.write('%d\t%d\t%d\t%d\t%.2f\n' % (q, qv_tot[q], e, qv_mm[q], emp))
    print('records %d aligned_query_bases %d' % (nrec, nbases), file=sys.stderr)


if __name__ == '__main__':
    main()
