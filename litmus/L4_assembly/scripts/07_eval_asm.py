#!/usr/bin/env python3
"""Evaluate a diploid (hap1/hap2) assembly of one chromosome against trio-phased truth haplotypes.
Metrics per assembly haplotype file and combined:
  n_ctg, total_bp, NG50 (G = mean truth haplotype length), NGA50 (primary minimap2 asm5 alignment blocks >=1 kb of
  each contig against the truth haplotype (pat or mat) to which it aligns best), cov999_{pat,mat}: fraction of each truth haplotype covered by an
  assembly alignment whose identity in every 10-kb target window is >= 99.9% (non-competitive: all contigs
  aligned to each truth haplotype separately; identity = matches / (matches + mismatches + indel bases)),
  QV (yak qv, Illumina k-mers of the same chromosome), switch/hamming error (yak trioeval, parental Illumina k-mers).
Usage: 07_eval_asm.py --asm DIR --pat pat.fa --mat mat.fa --child child.yak --patyak p.yak --matyak m.yak
                      --label NAME --out metrics.tsv [--threads 3] [--files hap1.fa,hap2.fa]
"""
import argparse, subprocess, os, re, tempfile, shutil, sys, json
import numpy as np
ap = argparse.ArgumentParser()
for k in ['asm', 'pat', 'mat', 'child', 'patyak', 'matyak', 'label', 'out']: ap.add_argument('--' + k, required=True)
ap.add_argument('--threads', default='3'); ap.add_argument('--files', default='hap1.fa,hap2.fa')
ap.add_argument('--win', type=int, default=10000); ap.add_argument('--keep-paf', default='')
a = ap.parse_args()
tmp = tempfile.mkdtemp(prefix='scratch.', dir=os.environ.get('TMPDIR', os.getcwd()))
def run(cmd, out=None):
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    if out: open(out, 'w').write(r.stdout)
    return r.stdout
def fa_lens(fn):
    L = {}; name = None
    for line in open(fn):
        if line.startswith('>'): name = line[1:].split()[0]; L[name] = 0
        elif name: L[name] += len(line.strip())
    return L
def n50(lengths, G):
    s = 0
    for x in sorted(lengths, reverse=True):
        s += x
        if s >= G / 2: return x
    return 0
cs_re = re.compile(r'(:\d+|\*[a-z][a-z]|[+\-][a-z]+)')
def cov999(paf_fn, tlen, win):
    covered = np.zeros(tlen, dtype=bool); aligned = np.zeros(tlen, dtype=bool)
    for line in open(paf_fn):
        f = line.rstrip('\n').split('\t')
        if 'tp:A:P' not in f[12:]: continue
        if int(f[10]) < 1000: continue
        cs = next((x[5:] for x in f[12:] if x.startswith('cs:Z:')), None)
        if cs is None: continue
        t = int(f[7]); tend = int(f[8])
        aligned[t:tend] = True
        # per-window accumulation for this alignment
        nwin = (tend - 1) // win - t // win + 1; w0 = t // win
        m = np.zeros(nwin); e = np.zeros(nwin); tc = np.zeros(nwin)
        pos = t
        for op in cs_re.findall(cs):
            c = op[0]
            if c == ':':
                n = int(op[1:]); p = pos
                while n > 0:
                    wi = p // win - w0; room = (p // win + 1) * win - p; k = min(room, n)
                    m[wi] += k; tc[wi] += k; p += k; n -= k
                pos = p
            elif c == '*':
                wi = pos // win - w0; e[wi] += 1; tc[wi] += 1; pos += 1
            elif c == '+':
                wi = min(pos // win - w0, nwin - 1); e[wi] += len(op) - 1
            elif c == '-':
                n = len(op) - 1
                for _ in range(n):
                    wi = pos // win - w0; e[wi] += 1; tc[wi] += 1; pos += 1
        ok = (m / np.maximum(m + e, 1)) >= 0.999
        for i in np.nonzero(ok & (tc > 0))[0]:
            ws = max(t, (w0 + i) * win); we = min(tend, (w0 + i + 1) * win)
            covered[ws:we] = True
    return covered.mean(), aligned.mean()

def qwin_errors(paf_fn, qwin=20000):
    """per (contig, query window): [aligned_query_bases, errors] from primary alignments (cs walk in query coords)"""
    out = {}
    for line in open(paf_fn):
        f = line.rstrip('\n').split('\t')
        if 'tp:A:P' not in f[12:] or int(f[10]) < 1000: continue
        cs = next((x[5:] for x in f[12:] if x.startswith('cs:Z:')), None)
        if cs is None: continue
        q, qs, qe, strand = f[0], int(f[2]), int(f[3]), f[4]
        d = out.setdefault(q, {}); off = 0
        def add(n, err):
            # query interval [off, off+n) in alignment order -> contig coordinates, split at window boundaries
            nonlocal off
            a0 = qs + off if strand == '+' else qe - off - n
            a1 = a0 + n; pos = a0; first = True
            while pos < a1:
                w = pos // qwin; nxt = min(a1, (w + 1) * qwin)
                v = d.setdefault(w, [0, 0]); v[0] += nxt - pos
                if first and err: v[1] += err; first = False
                pos = nxt
            off += n
        for op in cs_re.findall(cs):
            c = op[0]
            if c == ':': add(int(op[1:]), 0)
            elif c == '*': add(1, 1)
            elif c == '+': add(len(op) - 1, len(op) - 1)
            elif c == '-':
                w = (qs + off) // qwin if strand == '+' else max(0, qe - off - 1) // qwin
                v = d.setdefault(w, [0, 0]); v[1] += len(op) - 1
    return out
def truth_phasing(qe_pat, qe_mat, contigs, min_diff=3, min_aln=0.5, qwin=20000):
    """hamming = minority informative windows / informative windows; switch = haplotype changes between consecutive
    informative windows / (informative windows - contigs with >=1 informative window)"""
    inf = mino = sw = adj = 0
    for q in contigs:
        calls = []
        for w in sorted(set(qe_pat.get(q, {})) | set(qe_mat.get(q, {}))):
            p = qe_pat.get(q, {}).get(w); m = qe_mat.get(q, {}).get(w)
            if not p or not m or p[0] < min_aln * qwin or m[0] < min_aln * qwin: continue
            if p[1] + min_diff <= m[1]: calls.append('P')
            elif m[1] + min_diff <= p[1]: calls.append('M')
        if not calls: continue
        inf += len(calls); mino += min(calls.count('P'), calls.count('M'))
        sw += sum(1 for i in range(1, len(calls)) if calls[i] != calls[i - 1]); adj += len(calls) - 1
    return (mino / inf if inf else float('nan')), (sw / adj if adj else float('nan')), inf

def yak_qv(fa):
    out = run(['yak', 'qv', '-t', a.threads, '-p', a.child, fa])
    q = [l.split('\t') for l in out.split('\n') if l.startswith('QV')]
    return float(q[-1][2]) if q else float('nan'), float(q[-1][1]) if q else float('nan')
def yak_trio(fa):
    out = run(['yak', 'trioeval', '-t', a.threads, a.patyak, a.matyak, fa])
    res = {}
    for l in out.split('\n'):
        f = l.split('\t')
        if f[0] == 'W': res['switch_err'] = float(f[3]); res['switch_n'] = int(f[2])
        if f[0] == 'H': res['hamming_err'] = float(f[3]); res['hamming_n'] = int(f[2])
    return res
try:
    tl = {'pat': sum(fa_lens(a.pat).values()), 'mat': sum(fa_lens(a.mat).values())}
    G = (tl['pat'] + tl['mat']) / 2
    files = [os.path.join(a.asm, x) for x in a.files.split(',')]
    allfa = os.path.join(tmp, 'all.fa'); open(allfa, 'w').write(''.join(open(x).read() for x in files))
    rows = []
    # non-competitive completeness at >= 99.9% identity (all contigs vs each truth haplotype)
    comp = {}; blk = {'pat': {}, 'mat': {}}; score = {'pat': {}, 'mat': {}}
    for h in ('pat', 'mat'):
        paf = os.path.join(tmp, f'{h}.paf')
        run(['minimap2', '-cx', 'asm5', '--cs', '-t', a.threads, getattr(a, h), allfa], out=paf)
        if a.keep_paf: shutil.copy(paf, f'{a.keep_paf}.{h}.paf')
        comp[h] = cov999(paf, tl[h], a.win)
        for line in open(paf):
            f = line.split('\t')
            if 'tp:A:P' not in line or int(f[10]) < 1000: continue
            blk[h].setdefault(f[0], []).append(int(f[3]) - int(f[2])); score[h][f[0]] = score[h].get(f[0], 0) + int(f[9])
    qe = {h: qwin_errors(os.path.join(tmp, f'{h}.paf')) for h in ('pat', 'mat')}
    best = {}
    for q in set(score['pat']) | set(score['mat']):
        hb = 'pat' if score['pat'].get(q, 0) >= score['mat'].get(q, 0) else 'mat'; best[q] = blk[hb].get(q, [])
    for fn in files + [allfa]:
        tag = os.path.basename(fn).replace('.fa', '') if fn != allfa else 'hap1+hap2'
        L = fa_lens(fn); g = G if fn != allfa else 2 * G
        r = dict(label=a.label, file=tag, n_ctg=len(L), total_bp=sum(L.values()), NG50=n50(L.values(), g))
        blocks = [b for q in L for b in best.get(q, [])]
        r['NGA50'] = n50(blocks, g); r['n_blocks_ge10kb'] = sum(1 for b in blocks if b >= 10000)
        if fn != allfa:
            r['QV'], r['QV_raw'] = yak_qv(fn); r.update(yak_trio(fn))
            r['hamming_truth'], r['switch_truth'], r['phase_windows'] = truth_phasing(qe['pat'], qe['mat'], list(L))
        else:
            r['cov999_pat'], r['aln_pat'] = comp['pat']; r['cov999_mat'], r['aln_mat'] = comp['mat']
            r['cov999_mean'] = (comp['pat'][0] + comp['mat'][0]) / 2
        rows.append(r)
    keys = ['label', 'file', 'n_ctg', 'total_bp', 'NG50', 'NGA50', 'n_blocks_ge10kb', 'QV', 'QV_raw', 'switch_err', 'switch_n',
            'hamming_err', 'hamming_n', 'hamming_truth', 'switch_truth', 'phase_windows', 'cov999_pat', 'cov999_mat', 'cov999_mean', 'aln_pat', 'aln_mat']
    new = not os.path.exists(a.out)
    with open(a.out, 'a') as f:
        if new: f.write('\t'.join(keys) + '\n')
        for r in rows: f.write('\t'.join(str(round(r[k], 6) if isinstance(r.get(k), float) else r.get(k, '')) for k in keys) + '\n')
    for r in rows: print(json.dumps(r))
finally:
    shutil.rmtree(tmp, ignore_errors=True)
