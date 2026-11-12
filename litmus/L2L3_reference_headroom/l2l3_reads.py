#!/usr/bin/env python3
"""L2/L3 pilot: per-read alignment metrics for one reference (haploid) or the per-read
best of two haplotype references (diploid 'oracle' / personalized pair).

Input BAMs must be minimap2 output in INPUT ORDER (unsorted; identical FASTQ input for
both haplotypes), so that records of each read are consecutive and the two BAMs can be
walked in lockstep.

  haploid : l2l3_reads.py --bam X.bam --out X.reads.npz
  diploid : l2l3_reads.py --bam H1.bam --bam2 H2.bam --out D.reads.npz \
                [--chosen1 D.h1.bam --chosen2 D.h2.bam]

Per read (each mate separately for paired-end) we record: read-name CRC, mate, read
length, mapped, MAPQ, contig id, start, end, NM, alignment columns (M+I+D), aligned
query bases, soft clips (left/right), AS, proper-pair flag, number of supplementary
records and their aligned query bases, minimap2 'de' tag, chosen haplotype, number of gap
(I/D) events and bases in indels >= 50 bp (so that SV-spanning reads can be separated).
Best-of-two rule: per template (read or read pair) the haplotype with the larger sum of
AS over all primary+supplementary records; ties -> larger summed MAPQ -> hap1.
"""
import argparse, sys, zlib
import numpy as np
import pysam

FIELDS = [('crc', np.uint32), ('mate', np.int8), ('qlen', np.int32), ('mapped', np.int8),
          ('mapq', np.uint8), ('tid', np.int32), ('pos', np.int32), ('end', np.int32),
          ('nm', np.int32), ('cols', np.int32), ('qaln', np.int32), ('clip_l', np.int32),
          ('clip_r', np.int32), ('as_', np.int32), ('proper', np.int8), ('nsupp', np.int16),
          ('supp_qaln', np.int32), ('de', np.float32), ('hap', np.int8), ('mate_tid', np.int32),
          ('gap_ev', np.int32), ('big_gap_bases', np.int32)]


def groups(bam):
    """yield lists of records belonging to one template (consecutive same qname)"""
    cur, name = [], None
    for r in bam.fetch(until_eof=True):
        if r.is_secondary:
            continue
        if r.query_name != name and cur:
            yield cur
            cur = []
        name = r.query_name
        cur.append(r)
    if cur:
        yield cur


def cig_stats(r, big=50):
    """M/=/X, I, D, left clip, right clip, gap events, bases in indels >= big"""
    m = i = d = 0
    cl = cr = 0
    gev = bigb = 0
    ct = r.cigartuples or []
    for k, (op, n) in enumerate(ct):
        if op in (0, 7, 8):
            m += n
        elif op == 1:
            i += n
            gev += 1
            if n >= big:
                bigb += n
        elif op == 2:
            d += n
            gev += 1
            if n >= big:
                bigb += n
        elif op in (4, 5):
            if k == 0:
                cl = n
            else:
                cr = n
    return m, i, d, cl, cr, gev, bigb


def template_records(grp):
    """split a template group into per-mate (primary, [supplementaries])"""
    out = {}
    for r in grp:
        mate = 1 if r.is_read1 else (2 if r.is_read2 else 0)
        e = out.setdefault(mate, [None, []])
        if r.is_supplementary:
            e[1].append(r)
        else:
            e[0] = r
    return out


def score(grp):
    s = mq = 0
    for r in grp:
        if not r.is_unmapped:
            s += r.get_tag('AS') if r.has_tag('AS') else 0
            mq += r.mapping_quality
    return s, mq


def read_rows(grp, hap):
    rows = []
    for mate, (p, sup) in sorted(template_records(grp).items()):
        if p is None:
            continue
        qlen = p.infer_read_length() or (len(p.query_sequence) if p.query_sequence else 0)
        crc = zlib.crc32((p.query_name + '/' + str(mate)).encode())
        if p.is_unmapped:
            rows.append((crc, mate, qlen, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, np.nan, hap, p.next_reference_id if p.is_paired else -1, 0, 0))
            continue
        m, i, d, cl, cr, gev, bigb = cig_stats(p)
        nm = p.get_tag('NM') if p.has_tag('NM') else 0
        as_ = p.get_tag('AS') if p.has_tag('AS') else 0
        de = p.get_tag('de') if p.has_tag('de') else np.nan
        sq = 0
        for s in sup:
            sm, si, _, _, _, _, _ = cig_stats(s)
            sq += sm + si
        rows.append((crc, mate, qlen, 1, p.mapping_quality, p.reference_id, p.reference_start,
                     p.reference_end, nm, m + i + d, m + i, cl, cr, as_,
                     1 if p.is_proper_pair else 0, len(sup), sq, de, hap,
                     p.next_reference_id if p.is_paired else -1, gev, bigb))
    return rows


def to_struct(rows):
    arr = np.zeros(len(rows), dtype=FIELDS)
    if rows:
        cols = list(zip(*rows))
        for (name, _), c in zip(FIELDS, cols):
            arr[name] = np.asarray(c)
    return arr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--bam2')
    ap.add_argument('--out', required=True)
    ap.add_argument('--chosen1', help='BAM of records of templates assigned to hap1')
    ap.add_argument('--chosen2', help='BAM of records of templates assigned to hap2')
    ap.add_argument('--prefix1', default='h1_')
    ap.add_argument('--prefix2', default='h2_')
    a = ap.parse_args()

    b1 = pysam.AlignmentFile(a.bam, 'rb', check_sq=False)
    rows = []
    if not a.bam2:
        for g in groups(b1):
            rows.extend(read_rows(g, 0))
        np.savez_compressed(a.out, reads=to_struct(rows),
                            refnames=np.array(b1.references), reflens=np.array(b1.lengths))
        print('reads', len(rows), file=sys.stderr)
        return
    b2 = pysam.AlignmentFile(a.bam2, 'rb', check_sq=False)
    n1 = len(b1.references)
    o1 = o2 = None
    if a.chosen1:
        o1 = pysam.AlignmentFile(a.chosen1, 'wb', template=b1)
        o2 = pysam.AlignmentFile(a.chosen2, 'wb', template=b2)
    nchoose = [0, 0, 0]
    for g1, g2 in zip(groups(b1), groups(b2)):
        if g1[0].query_name != g2[0].query_name:
            raise SystemExit('BAMs not in identical input order at %s / %s' % (g1[0].query_name, g2[0].query_name))
        s1, m1 = score(g1)
        s2, m2 = score(g2)
        if (s2, m2) > (s1, m1):
            hap, g, o = 2, g2, o2
        else:
            hap, g, o = 1, g1, o1
        nchoose[hap] += 1
        if s1 == s2:
            nchoose[0] += 1
        rows.extend(read_rows(g, hap))
        if o is not None:
            for r in g:
                if not r.is_unmapped:
                    o.write(r)
    if o1 is not None:
        o1.close(); o2.close()
    refn = np.array([a.prefix1 + n for n in b1.references] + [a.prefix2 + n for n in b2.references])
    refl = np.array(list(b1.lengths) + list(b2.lengths))
    arr = to_struct(rows)
    # tids of hap2 rows are shifted into the combined namespace
    m2 = (arr['hap'] == 2) & (arr['tid'] >= 0)
    arr['tid'][m2] += n1
    m2b = (arr['hap'] == 2) & (arr['mate_tid'] >= 0)
    arr['mate_tid'][m2b] += n1
    np.savez_compressed(a.out, reads=arr, refnames=refn, reflens=refl,
                        nchoose=np.array(nchoose))
    print('templates: ties %d hap1 %d hap2 %d' % tuple(nchoose), file=sys.stderr)


if __name__ == '__main__':
    main()
