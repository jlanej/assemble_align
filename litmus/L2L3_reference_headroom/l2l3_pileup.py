#!/usr/bin/env python3
"""L2/L3 pilot: count REF/ALT (GRCh38 sense) alleles at truth heterozygous SNVs in one
coordinate-sorted, indexed BAM.

  l2l3_pileup.py --bam X.bam --fasta X.fa --sites sites.tsv --minbq 20 --out counts.tsv.gz

sites.tsv (header): site_id contig pos0 ref_x alt_x
  ref_x/alt_x are the GRCh38 REF/ALT alleles expressed in the forward orientation of
  the BAM's reference X (complemented when X is reverse-oriented relative to GRCh38).
Uses samtools mpileup (-B, no BAQ; read-pair overlap removal on; secondary/QC-fail/dup
filtered); counts are reported at MAPQ>=0 and MAPQ>=20.
"""
import argparse, gzip, os, subprocess, sys, tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--sites', required=True)
    ap.add_argument('--minbq', type=int, default=20)
    ap.add_argument('--mq', type=int, default=20)
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    sites = {}
    tmp = tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False, dir=os.path.dirname(os.path.abspath(a.out)))
    with open(a.sites) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {k: i for i, k in enumerate(hdr)}
        for l in fh:
            f = l.rstrip('\n').split('\t')
            c, p = f[ix['contig']], int(f[ix['pos0']])
            sites.setdefault((c, p), []).append((f[ix['site_id']], f[ix['ref_x']].upper(), f[ix['alt_x']].upper()))
            tmp.write('%s\t%d\t%d\n' % (c, p, p + 1))
    tmp.close()
    cmd = [a.samtools, 'mpileup', '-B', '-Q', str(a.minbq), '-q', '0', '-s', '-d', '0',
           '--ff', 'UNMAP,SECONDARY,QCFAIL,DUP', '--no-output-ins', '--no-output-ins',
           '--no-output-del', '--no-output-del', '--no-output-ends',
           '-l', tmp.name, '-f', a.fasta, a.bam]
    ctgs = {k[0] for k in sites}
    if len(ctgs) == 1:
        cmd[-1:-1] = ['-r', next(iter(ctgs))]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    res = {}
    for l in p.stdout:
        f = l.rstrip('\n').split('\t')
        c, pos1, xref = f[0], int(f[1]), f[2].upper()
        key = (c, pos1 - 1)
        if key not in sites:
            continue
        bases, mqs = f[4], f[6]
        if len(bases) != len(mqs):
            continue
        for sid, rx, ax in sites[key]:
            cnt = [0] * 8  # ref0 alt0 oth0 del0 ref20 alt20 oth20 del20
            for b, q in zip(bases, mqs):
                if b in '.,':
                    b = xref
                elif b in '*#':
                    b = '*'
                else:
                    b = b.upper()
                k = 0 if b == rx else (1 if b == ax else (3 if b == '*' else 2))
                cnt[k] += 1
                if ord(q) - 33 >= a.mq:
                    cnt[4 + k] += 1
            res[sid] = cnt + ['R' if xref == rx else ('A' if xref == ax else 'O')]
    p.wait()
    os.unlink(tmp.name)
    with gzip.open(a.out, 'wt') as fo:
        fo.write('site_id\tref_q0\talt_q0\toth_q0\tdel_q0\tref_q%d\talt_q%d\toth_q%d\tdel_q%d\txref_class\n' % ((a.mq,) * 4))
        for v in sites.values():
            for sid, _, _ in v:
                fo.write(sid + '\t' + '\t'.join(map(str, res.get(sid, [0] * 8 + ['.']))) + '\n')
    print('sites %d with_coverage %d' % (sum(len(v) for v in sites.values()), len(res)), file=sys.stderr)


if __name__ == '__main__':
    main()
