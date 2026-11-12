#!/usr/bin/env python3
"""Score assemblies and assembly-QC outputs against the HG01943 and HG01103 flagged-case regression tests.

Subcommands
  asm-hg01943  --pat P.fa[,P2.fa] --mat M.fa[,M2.fa] --junction-kmers hg01943_junction_31mers.tsv --ref-window win.fa
               [--hapmers-mat mat.meryl --hapmers-pat pat.meryl] [--threads 2] --out result.json
      win.fa = GRCh38 chr2:46,400,001-48,650,000 as one sequence (samtools faidx <GRCh38> chr2:46400001-48650000).
  asm-hg01103  --hap1 H1.fa[.gz] --hap2 H2.fa[.gz] --junction-kmers hg01103_junction_31mers.tsv --out result.json
  qc-hg01943   --bed calls.bed[,calls2.bed] --out result.json
      calls in HPRC release-2 HG01943 contig coordinates (PanSN or bare accession; 0-based BED); 4th column = label.
Requires: python>=3.9, pyahocorasick; minimap2 on PATH for asm-hg01943; meryl-lookup only for the optional hapmer test.
"""
import argparse, gzip, io, json, os, re, shutil, subprocess, sys, tempfile, collections
SEG = (46970725, 48057041); WIN0 = 46400000; LWIN = (46850001, 48180000)
def open_any(p):
    if p.startswith(('http://', 'https://')):
        import urllib.request
        fh = urllib.request.urlopen(p, timeout=600)
    else:
        fh = open(p, 'rb')
    return io.TextIOWrapper(gzip.GzipFile(fileobj=fh) if p.endswith('.gz') else fh, encoding='ascii', errors='replace')
def read_fasta(paths):
    for p in paths:
        name, buf = None, []
        for line in open_any(p):
            if line.startswith('>'):
                if name is not None: yield p, name, ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
        if name is not None: yield p, name, ''.join(buf).upper()
def kmer_automaton(tsv):
    import ahocorasick
    comp = str.maketrans('ACGTN', 'TGCAN'); A = ahocorasick.Automaton(); groups = collections.Counter()
    for line in open(tsv):
        n, s = line.split()[:2]; s = s.upper(); g = n.split('_')[0]; groups[g] += 1
        A.add_word(s, (n, g, len(s))); A.add_word(s.translate(comp)[::-1], (n, g, len(s)))
    A.make_automaton(); return A, groups
def junction_scan(paths, tsv):
    A, groups = kmer_automaton(tsv); found = collections.defaultdict(set); where = collections.defaultdict(collections.Counter); lens = {}
    for p, name, seq in read_fasta(paths):
        lens[name] = len(seq)
        for end, (n, g, k) in A.iter(seq):
            found[g].add(n); where[g][(name, end)] += 1
    ctx = {}
    for g, c in where.items():
        cnt = collections.Counter(nm for nm, _ in c)
        nm = cnt.most_common(1)[0][0]; pos = sorted(e for (n2, e) in c if n2 == nm); p = pos[len(pos) // 2]; L = lens[nm]
        ctx[g] = dict(contig=nm, contig_length=L, distance_to_contig_end=min(p, L - p),
                      context='short contig (<100 kb)' if L < 100000 else ('contig end (<100 kb)' if min(p, L - p) < 100000 else 'internal'))
    return {g: len(found[g]) for g in groups}, groups, ctx
def cmd_asm_hg01943(a):
    pat, mat = a.pat.split(','), a.mat.split(',')
    res = {'test': 'R1 assembler: HG01943 chr2 1.09-Mb maternal tandem duplication'}
    jp, groups, _ = junction_scan(pat, a.junction_kmers); jm, _, ctxm = junction_scan(mat, a.junction_kmers)
    nJ = groups['J']; res['junction_kmers_total'] = nJ
    res['junction_kmers_in_maternal'] = jm.get('J', 0); res['junction_kmers_in_paternal'] = jp.get('J', 0)
    res['junction_context_maternal'] = ctxm.get('J')
    tmp = tempfile.mkdtemp(prefix='score_'); cn = {}
    try:
        for lab, paths in [('maternal', mat), ('paternal', pat)]:
            fa = os.path.join(tmp, lab + '.fa')
            with open(fa, 'w') as f:   # contig-level evaluation: split scaffolds at runs of >=10 N
                for p, name, seq in read_fasta(paths):
                    for k, piece in enumerate(re.split('N{10,}', seq)):
                        if len(piece) >= 1000: f.write(f'>{name}__piece{k}\n{piece}\n')
            paf = subprocess.run(['minimap2', '-t', str(a.threads), '-cx', 'asm5', '--secondary=no', a.ref_window, fa], capture_output=True, text=True, check=True).stdout
            nb = (2250000 + 9999) // 10000; dep = [0] * nb; gaps = 0; names = set()
            for line in paf.splitlines():
                t = line.split('\t')
                if int(t[10]) < 10000: continue
                names.add(t[0].split('__piece')[0]); s, e = int(t[7]), int(t[8])
                for b in range(s // 10000, (e - 1) // 10000 + 1): dep[b] += min(e, (b + 1) * 10000) - max(s, b * 10000)
            for p, name, seq in read_fasta(paths):
                if name in names: gaps += len(re.findall('N{10,}', seq))
            good = tot = 0; disc = []
            for b in range(nb):
                s1, e1 = WIN0 + b * 10000 + 1, WIN0 + (b + 1) * 10000
                if s1 < LWIN[0] or e1 > LWIN[1] or (s1 <= SEG[0] <= e1) or (s1 <= SEG[1] <= e1): continue
                exp = 2 if (lab == 'maternal' and SEG[0] < s1 and e1 < SEG[1]) else 1
                obs = round(dep[b] / 10000); tot += 1; good += obs == exp
                if obs != exp: disc.append((s1, e1, exp, obs))
            cn[lab] = dict(bins=tot, concordant=good, fraction=round(good / max(1, tot), 3), gaps_in_aligned_contigs=gaps, discordant_bins_first10=disc[:10])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    res['copy_number_vs_expected'] = cn
    if (a.hapmers_mat or a.hapmers_pat) and not shutil.which('meryl-lookup'):
        sys.exit('hapmer test requested but meryl-lookup is not on PATH')
    if a.hapmers_mat and a.hapmers_pat:
        bad = []
        tmp2 = tempfile.mkdtemp(prefix='score_hm_')
        for lab, paths in [('paternal', pat), ('maternal', mat)]:
            pf = os.path.join(tmp2, lab + '.pieces.fa')
            with open(pf, 'w') as f:   # scaffold pieces between N-runs are scored separately
                for _, name, seq in read_fasta(paths):
                    for k, piece in enumerate(re.split('N{10,}', seq)):
                        if len(piece) >= 1000: f.write(f'>{name}__piece{k}\n{piece}\n')
            for p in [pf]:
                out = subprocess.run(['meryl-lookup', '-existence', '-sequence', p, '-mers', a.hapmers_mat, a.hapmers_pat], capture_output=True, text=True).stdout
                for line in out.splitlines():
                    t = line.split('\t'); m, pt = int(t[3]), int(t[5])
                    wrong = (m > 3 * pt and m >= 100) if lab == 'paternal' else (pt > 3 * m and pt >= 100)
                    if wrong: bad.append(dict(haplotype=lab, contig=t[0], maternal_hapmers=m, paternal_hapmers=pt))
        shutil.rmtree(tmp2, ignore_errors=True)
        res['trio_inconsistent_contigs'] = bad
    ok_j = res['junction_kmers_in_maternal'] >= 0.9 * nJ and res['junction_kmers_in_paternal'] == 0
    ok_cn = all(v['fraction'] >= 0.95 for v in cn.values()); ok_gap = all(v['gaps_in_aligned_contigs'] == 0 for v in cn.values())
    ok_trio = not res.get('trio_inconsistent_contigs')
    res['criteria'] = dict(junction_on_maternal_only=ok_j, copy_number_concordance_ge_0_95=ok_cn, no_gap_joins=ok_gap, trio_consistent=ok_trio)
    res['verdict'] = 'PASS' if all(res['criteria'].values()) else 'FAIL'
    json.dump(res, open(a.out, 'w'), indent=1); print(json.dumps({k: v for k, v in res.items() if k != 'copy_number_vs_expected'}, indent=1)); print(json.dumps(cn, indent=1))
def cmd_asm_hg01103(a):
    res = {'test': 'R3 assembler: HG01103 sub-clonal chr2/chr14 chromothripsis (16 junctions)'}
    per = {}
    for lab, p in [('hap1', a.hap1), ('hap2', a.hap2)]:
        f, groups, ctx = junction_scan(p.split(','), a.junction_kmers)
        per[lab] = {g: dict(kmers_found=f.get(g, 0), kmers=groups[g], present=f.get(g, 0) >= 0.8 * groups[g], **(ctx.get(g) or {})) for g in sorted(groups, key=lambda x: int(x[1:]))}
    res['per_junction'] = per
    n1 = sum(v['present'] for v in per['hap1'].values()); n2 = sum(v['present'] for v in per['hap2'].values())
    one = max(n1, n2); which = 'hap1' if n1 >= n2 else 'hap2'
    internal = sum(1 for v in per[which].values() if v['present'] and v.get('context') == 'internal')
    res['summary'] = dict(junctions_hap1=n1, junctions_hap2=n2, junctions_any=sum(per['hap1'][g]['present'] or per['hap2'][g]['present'] for g in per['hap1']),
                          derivative_haplotype=which, junctions_on_one_haplotype=one, junctions_internal_to_contigs=internal)
    res['criteria'] = dict(all_16_on_one_haplotype=(one == 16 and min(n1, n2) == 0), threaded_ge_15_internal=(internal >= 15))
    res['verdict'] = 'PASS' if all(res['criteria'].values()) else 'FAIL'
    json.dump(res, open(a.out, 'w'), indent=1); print(json.dumps({k: v for k, v in res.items() if k != 'per_junction'}, indent=1))
# QC targets for HG01943 in HPRC r2 contig coordinates (1-based, inclusive)
P, M = 'CM089354.1', 'CM089362.1'
TARGETS = {   # target -> list of sub-targets; each sub-target = list of alternative intervals (contig, start, end)
 'T1_unsupported_joins': [[(P, 47394431, 47398531)], [(P, 47763590, 47767690)], [(P, 48219141, 48223241)], [(P, 48520660, 48524760)]],
 'T2_trio_inconsistent_piece_B': [[(P, 47396532, 47765590)]],
 'T3_missing_copy_collapse': [[(M, 47368148, 47454263)], [(P, 48102850, 48201244), (M, 47684402, 47782807)]],
 'T4_junction_absent_M_segment_ends': [[(M, 46992979, 47002979)], [(M, 48080387, 48090387)]],
 'T5_false_duplications': [[(P, 47375929, 47396431), (P, 47765691, 47786203)], [(P, 48202383, 48221142), (P, 48522761, 48541520)]],
}
RULE = {'T1_unsupported_joins': (1e-9, 'all'), 'T2_trio_inconsistent_piece_B': (0.5, 'all'), 'T3_missing_copy_collapse': ('5kb', 'all'),
        'T4_junction_absent_M_segment_ends': (1e-9, 'any'), 'T5_false_duplications': (0.25, 'all')}
LOCUS = {P: (46450001, 49420000), M: (46420001, 48680000), 'JBHIIL010000075.1': (1, 230135)}
GAPLAB = {'NNN', 'GAP', 'N'}
def cmd_qc_hg01943(a):
    calls = []
    for p in a.bed.split(','):
        for line in open(p):
            if line.startswith(('#', 'track', 'browser')) or not line.strip(): continue
            t = line.rstrip('\n').split('\t'); c = t[0].split('#')[-1]
            calls.append((c, int(t[1]) + 1, int(t[2]), t[3] if len(t) > 3 else 'call'))
    calls = [x for x in calls if x[0] in LOCUS and x[1] <= LOCUS[x[0]][1] and x[2] >= LOCUS[x[0]][0]]
    def ov(c, s, e, x): return max(0, min(e, x[2]) - max(s, x[1]) + 1) if c == x[0] else 0
    hit = {}
    for T, subs in TARGETS.items():
        need, mode = RULE[T]; flags = []
        for alts in subs:
            ok = False
            for (c, s0, e0) in alts:
                cov = sum(ov(c, s0, e0, x) for x in calls if not (T == 'T1_unsupported_joins' and x[3].upper() in GAPLAB))
                ok = ok or (cov >= 5000 if need == '5kb' else cov / (e0 - s0 + 1) >= need)
            flags.append(ok)
        hit[T] = dict(subtargets_hit=sum(flags), subtargets=len(flags), hit=(all(flags) if mode == 'all' else any(flags)))
    def near(x):
        return any(x[0] == c and x[1] <= e + 5000 and x[2] >= s - 5000 for subs in TARGETS.values() for alts in subs for (c, s, e) in alts)
    off = [x for x in calls if x[3].upper() not in GAPLAB and not near(x)]
    res = dict(test='R2 assembly QC: HG01943 locus', n_calls_in_locus=len(calls), targets=hit, targets_hit=sum(v['hit'] for v in hit.values()),
               off_target_calls=len(off), off_target_bp=sum(x[2] - x[1] + 1 for x in off), off_target_labels=dict(collections.Counter(x[3] for x in off)))
    res['verdict'] = 'PASS' if res['targets_hit'] == 5 else 'FAIL'
    json.dump(res, open(a.out, 'w'), indent=1); print(json.dumps(res, indent=1))
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); sp = ap.add_subparsers(dest='cmd', required=True)
    s1 = sp.add_parser('asm-hg01943'); s1.add_argument('--pat', required=True); s1.add_argument('--mat', required=True); s1.add_argument('--junction-kmers', required=True)
    s1.add_argument('--ref-window', required=True); s1.add_argument('--hapmers-mat'); s1.add_argument('--hapmers-pat'); s1.add_argument('--threads', type=int, default=2); s1.add_argument('--out', required=True)
    s2 = sp.add_parser('asm-hg01103'); s2.add_argument('--hap1', required=True); s2.add_argument('--hap2', required=True); s2.add_argument('--junction-kmers', required=True); s2.add_argument('--out', required=True)
    s3 = sp.add_parser('qc-hg01943'); s3.add_argument('--bed', required=True); s3.add_argument('--out', required=True)
    a = ap.parse_args(); {'asm-hg01943': cmd_asm_hg01943, 'asm-hg01103': cmd_asm_hg01103, 'qc-hg01943': cmd_qc_hg01943}[a.cmd](a)
if __name__ == '__main__': main()
