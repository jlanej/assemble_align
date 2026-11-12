#!/usr/bin/env python3
"""L2/L3 reference-headroom pilot: assembly-based truth set and coordinate projection.

Builds, for one target chromosome (default chr20 of GRCh38), the NA19240 truth set
from minimap2 asm5 alignments (PAF with --cs) of the two phased haplotypes:
  * per-haplotype target->query projection arrays (dipcall-like: only target bases
    covered by exactly one >=MINLEN alignment of each haplotype are 'confident')
  * heterozygous / homozygous-ALT SNVs relative to GRCh38
  * for every het SNV, its coordinate + orientation on hap1, hap2 and CHM13
    (CHM13 via hap1->CHM13 and hap2->CHM13 alignments; kept only if both agree)

Usage:
  l2l3_truth.py --ref GRCh38.fa --chrom chr20 \
      --h1 hap1_vs_GRCh38.paf.gz --h2 hap2_vs_GRCh38.paf.gz \
      --h1c hap1_vs_CHM13.paf.gz --h2c hap2_vs_CHM13.paf.gz --chm13 CHM13.fa \
      --out truth_prefix
Outputs: <out>.het_snv.tsv.gz, <out>.homalt_snv.tsv.gz, <out>.confident.bed.gz,
         <out>.summary.json
"""
import argparse, gzip, json, re, sys
import numpy as np
import pysam

CS_RE = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def iter_paf(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for l in fh:
            f = l.rstrip('\n').split('\t')
            r = dict(qname=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                     tname=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8]),
                     nmatch=int(f[9]), blen=int(f[10]), mapq=int(f[11]))
            tags = {}
            for t in f[12:]:
                k, ty, v = t.split(':', 2)
                tags[k] = v
            r['tags'] = tags
            yield r


class TargetMap:
    """target-indexed projection of one haplotype onto one target chromosome"""

    def __init__(self, tlen):
        self.tlen = tlen
        self.diff = np.zeros(tlen + 1, np.int32)
        self.qid = np.full(tlen, -1, np.int16)
        self.qpos = np.full(tlen, -1, np.int32)
        self.strand = np.zeros(tlen, np.int8)
        self.snv = {}            # tpos -> alt base (target orientation)
        self.indel = []          # (tpos, type, len)
        self.qnames = []

    def add(self, r):
        if r['qname'] not in self.qnames:
            self.qnames.append(r['qname'])
        qid = self.qnames.index(r['qname'])
        fwd = r['strand'] == '+'
        t = r['ts']
        q = r['qs'] if fwd else r['qe'] - 1
        sgn = 1 if fwd else -1
        self.diff[r['ts']] += 1
        self.diff[r['te']] -= 1
        for op in CS_RE.findall(r['tags']['cs']):
            c = op[0]
            if c == ':':
                n = int(op[1:])
                self.qpos[t:t + n] = np.arange(q, q + sgn * n, sgn, dtype=np.int64)[:n]
                self.qid[t:t + n] = qid
                self.strand[t:t + n] = sgn
                t += n
                q += sgn * n
            elif c == '*':
                self.qpos[t] = q
                self.qid[t] = qid
                self.strand[t] = sgn
                self.snv[t] = op[2].upper()
                t += 1
                q += sgn
            elif c == '+':
                n = len(op) - 1
                self.indel.append((t, 'I', n))
                q += sgn * n
            elif c == '-':
                n = len(op) - 1
                self.indel.append((t, 'D', n))
                t += n
            elif c == '=':
                n = len(op) - 1
                self.qpos[t:t + n] = np.arange(q, q + sgn * n, sgn, dtype=np.int64)[:n]
                self.qid[t:t + n] = qid
                self.strand[t:t + n] = sgn
                t += n
                q += sgn * n
            else:
                raise ValueError('unsupported cs op ' + op)
        if t != r['te']:
            raise ValueError('cs walk mismatch %s %d %d' % (r['qname'], t, r['te']))

    def depth(self):
        return np.cumsum(self.diff)[:-1]


class QueryMap:
    """query-indexed projection (hap contig -> target genome) for selected contigs"""

    def __init__(self, qlens):
        self.q = {n: dict(tid=np.full(l, -1, np.int16), tpos=np.full(l, -1, np.int32),
                          strand=np.zeros(l, np.int8), diff=np.zeros(l + 1, np.int32))
                  for n, l in qlens.items()}
        self.tnames = []

    def add(self, r):
        if r['qname'] not in self.q:
            return
        d = self.q[r['qname']]
        if r['tname'] not in self.tnames:
            self.tnames.append(r['tname'])
        tid = self.tnames.index(r['tname'])
        fwd = r['strand'] == '+'
        sgn = 1 if fwd else -1
        t = r['ts']
        q = r['qs'] if fwd else r['qe'] - 1
        d['diff'][r['qs']] += 1
        d['diff'][r['qe']] -= 1
        for op in CS_RE.findall(r['tags']['cs']):
            c = op[0]
            if c in ':=':
                n = int(op[1:]) if c == ':' else len(op) - 1
                qq = np.arange(q, q + sgn * n, sgn, dtype=np.int64)[:n]
                d['tpos'][qq] = np.arange(t, t + n, dtype=np.int64)
                d['tid'][qq] = tid
                d['strand'][qq] = sgn
                t += n
                q += sgn * n
            elif c == '*':
                d['tpos'][q] = t
                d['tid'][q] = tid
                d['strand'][q] = sgn
                t += 1
                q += sgn
            elif c == '+':
                q += sgn * (len(op) - 1)
            elif c == '-':
                t += len(op) - 1

    def depth(self, name):
        return np.cumsum(self.q[name]['diff'])[:-1]


def load_target_map(paf, chrom, tlen, minlen):
    tm = TargetMap(tlen)
    n = 0
    for r in iter_paf(paf):
        if r['tname'] != chrom or r['blen'] < minlen or 'cs' not in r['tags']:
            continue
        if r['tags'].get('tp', 'P') not in ('P', 'I'):
            continue
        tm.add(r)
        n += 1
    return tm, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True)
    ap.add_argument('--chrom', default='chr20')
    ap.add_argument('--h1', required=True)
    ap.add_argument('--h2', required=True)
    ap.add_argument('--h1c')
    ap.add_argument('--h2c')
    ap.add_argument('--chm13')
    ap.add_argument('--minlen', type=int, default=10000)
    ap.add_argument('--indel-flank', type=int, default=5)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    fa = pysam.FastaFile(a.ref)
    seq = np.frombuffer(fa.fetch(a.chrom).upper().encode(), dtype='S1')
    tlen = len(seq)
    tm1, n1 = load_target_map(a.h1, a.chrom, tlen, a.minlen)
    tm2, n2 = load_target_map(a.h2, a.chrom, tlen, a.minlen)
    d1, d2 = tm1.depth(), tm2.depth()
    isacgt = np.isin(seq, [b'A', b'C', b'G', b'T'])
    conf = (d1 == 1) & (d2 == 1) & isacgt
    # near-indel mask
    near = np.zeros(tlen + 1, np.int32)
    fl = a.indel_flank
    for tm in (tm1, tm2):
        for t, ty, n in tm.indel:
            s = max(0, t - fl)
            e = min(tlen, t + (n if ty == 'D' else 0) + fl)
            near[s] += 1
            near[e] -= 1
    near = np.cumsum(near)[:-1] > 0

    cand = sorted(set(tm1.snv) | set(tm2.snv))
    rows_het, rows_hom = [], []
    for t in cand:
        if not conf[t] or tm1.qpos[t] < 0 or tm2.qpos[t] < 0:
            continue
        ref = seq[t].decode()
        a1 = tm1.snv.get(t, ref)
        a2 = tm2.snv.get(t, ref)
        if a1 == a2:
            if a1 != ref:
                rows_hom.append((t, ref, a1, int(near[t])))
            continue
        if a1 != ref and a2 != ref:
            continue  # het-alt (two ALT alleles) skipped
        alt = a1 if a1 != ref else a2
        rows_het.append((t, ref, alt, 1 if a1 == alt else 2,
                         int(tm1.qid[t]), int(tm1.qpos[t]), int(tm1.strand[t]),
                         int(tm2.qid[t]), int(tm2.qpos[t]), int(tm2.strand[t]), int(near[t])))

    # CHM13 projection through each haplotype
    chm = {}
    if a.h1c and a.h2c:
        for hap, paf, tm in ((1, a.h1c, tm1), (2, a.h2c, tm2)):
            # contig lengths from PAF records
            qlens = {}
            for r in iter_paf(paf):
                if r['qname'] in tm.qnames:
                    qlens[r['qname']] = r['qlen']
            qm = QueryMap(qlens)
            for r in iter_paf(paf):
                if r['blen'] < a.minlen or 'cs' not in r['tags']:
                    continue
                qm.add(r)
            dep = {n: qm.depth(n) for n in qlens}
            chm[hap] = (qm, dep, tm)
    out_het = []
    nchm_ok = 0
    for row in rows_het:
        t, ref, alt, althap, q1, p1, s1, q2, p2, s2, nr = row
        c_chr, c_pos, c_or = '.', -1, 0
        if chm:
            proj = []
            for hap, qid, qpos, s in ((1, q1, p1, s1), (2, q2, p2, s2)):
                qm, dep, tm = chm[hap]
                name = tm.qnames[qid]
                if name not in qm.q:
                    proj.append(None)
                    continue
                dd = qm.q[name]
                if dep[name][qpos] != 1 or dd['tpos'][qpos] < 0:
                    proj.append(None)
                    continue
                # orientation of GRCh38-forward relative to CHM13-forward = s(hap vs GRCh38) * s(hap vs CHM13)
                proj.append((qm.tnames[dd['tid'][qpos]], int(dd['tpos'][qpos]), int(s) * int(dd['strand'][qpos])))
            if proj[0] is not None and proj[0] == proj[1]:
                c_chr, c_pos, c_or = proj[0]
                nchm_ok += 1
        out_het.append((a.chrom, t, ref, alt, althap,
                        tm1.qnames[q1], p1, s1, tm2.qnames[q2], p2, s2, c_chr, c_pos, c_or, nr))

    with gzip.open(a.out + '.het_snv.tsv.gz', 'wt') as fo:
        fo.write('chrom\tpos0\tref\talt\talt_hap\th1_contig\th1_pos0\th1_orient\th2_contig\th2_pos0\th2_orient\t'
                 'chm13_chrom\tchm13_pos0\tchm13_orient\tnear_indel\n')
        for r in out_het:
            fo.write('\t'.join(map(str, r)) + '\n')
    with gzip.open(a.out + '.homalt_snv.tsv.gz', 'wt') as fo:
        fo.write('chrom\tpos0\tref\talt\tnear_indel\n')
        for t, ref, alt, nr in rows_hom:
            fo.write('%s\t%d\t%s\t%s\t%d\n' % (a.chrom, t, ref, alt, nr))
    # confident region BED
    with gzip.open(a.out + '.confident.bed.gz', 'wt') as fo:
        x = np.diff(np.concatenate([[0], conf.astype(np.int8), [0]]))
        st = np.where(x == 1)[0]
        en = np.where(x == -1)[0]
        for s, e in zip(st, en):
            fo.write('%s\t%d\t%d\n' % (a.chrom, s, e))
    summ = dict(chrom=a.chrom, tlen=int(tlen), h1_alignments=n1, h2_alignments=n2,
                h1_contigs=tm1.qnames, h2_contigs=tm2.qnames,
                confident_bp=int(conf.sum()), acgt_bp=int(isacgt.sum()),
                het_snv=len(out_het), het_snv_near_indel=int(sum(r[-1] for r in out_het)),
                homalt_snv=len(rows_hom), het_snv_with_chm13=nchm_ok,
                h1_indels=len(tm1.indel), h2_indels=len(tm2.indel))
    json.dump(summ, open(a.out + '.summary.json', 'w'), indent=1)
    print(json.dumps(summ))


if __name__ == '__main__':
    main()
