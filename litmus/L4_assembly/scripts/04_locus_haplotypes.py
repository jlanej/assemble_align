#!/usr/bin/env python3
"""Project GRCh38 locus intervals onto HPRC release-2 haplotypes with the per-assembly
minigraph-cactus chains (<asm>_vs_GRCh38.chain.gz; t = assembly contig, q = GRCh38) and fetch
the haplotype sequences with remote `samtools faidx` (bgzipped fa.gz + .fai + .gzi on S3).
Sequences are written in GRCh38 orientation (reverse-complemented when the best chain is on '-').

Usage: 04_locus_haplotypes.py <loci.tsv> <locus_assemblies.tsv> <outdir> [--truth SAMPLE,SAMPLE]
  loci.tsv: locus chrom core_start core_end flank genes win_start win_end core_len (0-based half-open)
  locus_assemblies.tsv: sample hap name fa_url chain_url
Outputs: <outdir>/projections.tsv, <outdir>/<locus>/pop.W.fa (panel), <outdir>/<locus>/truth.<sample>.<hap>.{W,C}.fa
"""
import sys, os, gzip, csv, subprocess, tempfile, shutil, urllib.request, argparse, collections

def parse_chains(path, chroms):
    chains = []
    with gzip.open(path, 'rt') as f:
        cur = None
        for line in f:
            if line.startswith('chain'):
                p = line.split()
                cur = dict(score=int(p[1]), t=p[2], tsize=int(p[3]), ts=int(p[5]), te=int(p[6]),
                           q=p[7], qsize=int(p[8]), qstrand=p[9], qs=int(p[10]), qe=int(p[11]), blocks=[])
                keep = cur['q'] in chroms
                tpos, qpos = cur['ts'], cur['qs']
                if keep: chains.append(cur)
                continue
            p = line.split()
            if not p or cur is None: continue
            if keep:
                size = int(p[0])
                # store block in forward GRCh38 coordinates
                if cur['qstrand'] == '+': fq0, fq1 = qpos, qpos + size
                else: fq0, fq1 = cur['qsize'] - qpos - size, cur['qsize'] - qpos
                cur['blocks'].append((fq0, fq1, tpos, size))
            if len(p) == 3:
                tpos += int(p[0]) + int(p[1]); qpos += int(p[0]) + int(p[2])
    return chains

def project(chains, chrom, s, e):
    """Return dict(t, tstart, tend, strand, aligned_frac, n_chains, ratio) for GRCh38 [s,e)."""
    L = e - s; per = []
    for c in chains:
        if c['q'] != chrom or c['qe'] <= 0: continue
        ov = 0; tmin = None; tmax = None; qmin = None; qmax = None
        for fq0, fq1, t0, size in c['blocks']:
            a, b = max(fq0, s), min(fq1, e)
            if a >= b: continue
            ov += b - a
            if c['qstrand'] == '+': ta, tb = t0 + (a - fq0), t0 + (b - fq0)
            else: ta, tb = t0 + (fq1 - b), t0 + (fq1 - a)
            tmin = ta if tmin is None else min(tmin, ta); tmax = tb if tmax is None else max(tmax, tb)
            qmin = a if qmin is None else min(qmin, a); qmax = b if qmax is None else max(qmax, b)
        if ov: per.append((ov, c, tmin, tmax, qmin, qmax))
    if not per: return None
    per.sort(key=lambda x: -x[0])
    ov, best, tmin, tmax, qmin, qmax = per[0]
    used = 1; tot = ov
    for ov2, c2, a2, b2, qa, qb in per[1:]:   # merge neighbouring chains on the same contig/strand
        if c2['t'] == best['t'] and c2['qstrand'] == best['qstrand'] and a2 > tmin - 2 * L and b2 < tmax + 2 * L:
            tmin, tmax = min(tmin, a2), max(tmax, b2); qmin, qmax = min(qmin, qa), max(qmax, qb); used += 1; tot += ov2
    # extend for unaligned ends of the GRCh38 interval (approximate, recorded)
    lext, rext = qmin - s, e - qmax
    if best['qstrand'] == '+': tmin -= lext; tmax += rext
    else: tmin -= rext; tmax += lext
    tmin = max(0, tmin); tmax = min(best['tsize'], tmax)
    return dict(t=best['t'], tstart=tmin, tend=tmax, strand=best['qstrand'], aligned_frac=round(tot / L, 4),
                n_chains=used, ratio=round((tmax - tmin) / L, 3), end_gap_left=lext, end_gap_right=rext)

def project_anchored(chains, chrom, s, e, alen=20000):
    """Project the unique flanking anchors [s,s+alen) and [e-alen,e) separately and return the haplotype
    interval between them (same contig and strand required); None if the anchors are inconsistent."""
    pl = project(chains, chrom, s, s + alen); pr = project(chains, chrom, e - alen, e)
    if not pl or not pr or pl['t'] != pr['t'] or pl['strand'] != pr['strand']: return None
    if pl['aligned_frac'] < 0.8 or pr['aligned_frac'] < 0.8: return None
    ts, te = (pl['tstart'], pr['tend']) if pl['strand'] == '+' else (pr['tstart'], pl['tend'])
    if te <= ts: return None
    return dict(t=pl['t'], tstart=ts, tend=te, strand=pl['strand'], aligned_frac=min(pl['aligned_frac'], pr['aligned_frac']),
                n_chains=-1, ratio=round((te - ts) / (e - s), 3), end_gap_left=pl['end_gap_left'], end_gap_right=pr['end_gap_right'])

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('loci'); ap.add_argument('asms'); ap.add_argument('outdir')
    ap.add_argument('--truth', default=''); ap.add_argument('--only', default=''); a = ap.parse_args()
    truth = set(x for x in a.truth.split(',') if x)
    loci = [l for l in csv.DictReader(open(a.loci), delimiter='\t') if not a.only or l['locus'] in a.only.split(',')]
    asms = list(csv.DictReader(open(a.asms), delimiter='\t'))
    chroms = set(l['chrom'] for l in loci)
    os.makedirs(a.outdir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='scratch.', dir=os.environ.get('TMPDIR', os.getcwd()))
    keys = ['sample', 'hap', 'name', 'locus', 'kind', 'grch38', 'status', 'contig', 'tstart', 'tend', 'strand',
            'aligned_frac', 'n_chains', 'ratio', 'end_gap_left', 'end_gap_right', 't']
    pt = os.path.join(a.outdir, 'projections.tsv')
    if not os.path.exists(pt): open(pt, 'w').write('\t'.join(keys) + '\n')
    def emit(r):
        with open(pt, 'a') as f: f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
    try:
        for asm in asms:
            cp = os.path.join(tmp, 'c.chain.gz'); urllib.request.urlretrieve(asm['chain_url'], cp)
            chains = parse_chains(cp, chroms)
            fai = urllib.request.urlopen(asm['fa_url'] + '.fai').read().decode().split('\n')
            names = {l.split('\t')[0].split('#')[-1]: l.split('\t')[0] for l in fai if l}
            clen = {l.split('\t')[0]: int(l.split('\t')[1]) for l in fai if l}
            is_truth = asm['sample'] in truth
            for l in loci:
                ivs = [('W', int(l['win_start']), int(l['win_end']))]
                if is_truth: ivs.append(('C', int(l['core_start']), int(l['core_end'])))
                for kind, s, e in ivs:
                    p = (project_anchored(chains, l['chrom'], s, e) if l.get('mode', 'overlap') == 'anchor'
                         else project(chains, l['chrom'], s, e))
                    row = dict(sample=asm['sample'], hap=asm['hap'], name=asm['name'], locus=l['locus'], kind=kind,
                               grch38=f"{l['chrom']}:{s}-{e}")
                    if p is None:
                        row.update(status='no_chain'); emit(row); continue
                    row.update(p); row['contig'] = names.get(p['t'], p['t'])
                    ok = 0.5 <= p['ratio'] <= 3.0 and (p['aligned_frac'] >= 0.3 or l.get('mode') == 'anchor')
                    row['status'] = 'ok' if ok else 'suspect'
                    # fetch
                    p['tend'] = min(p['tend'], clen.get(row['contig'], p['tend']))
                    reg = f"{row['contig']}:{p['tstart'] + 1}-{p['tend']}"
                    od = os.path.join(a.outdir, l['locus']); os.makedirs(od, exist_ok=True)
                    cmd = ['samtools', 'faidx'] + (['-i'] if p['strand'] == '-' else []) + [asm['fa_url'], reg]
                    seq = None
                    for attempt in range(3):
                        r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
                        if r.returncode == 0 and r.stdout.count('\n') > 1: seq = r.stdout; break
                    if seq is None:
                        row['status'] = 'fetch_failed'; emit(row); continue
                    hdr = f">{asm['sample']}#{asm['hap']}#{l['locus']}#{kind} {reg} strand={p['strand']} status={row['status']}\n"
                    body = ''.join(seq.split('\n')[1:])
                    fn = (f"truth.{asm['sample']}.{asm['hap']}.{kind}.fa" if is_truth else 'pop.W.fa')
                    with open(os.path.join(od, fn), 'a' if not is_truth else 'w') as fo:
                        fo.write(hdr + '\n'.join(body[i:i + 80] for i in range(0, len(body), 80)) + '\n')
                    emit(row)
            print(asm['name'], 'done', file=sys.stderr, flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == '__main__':
    main()
