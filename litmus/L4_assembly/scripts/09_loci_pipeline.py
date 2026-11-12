#!/usr/bin/env python3
"""Hard-locus recruitment -> hifiasm -> truth evaluation (one sample, HiFi reads).
Recruitment strategies (all drawn from the same candidate pool = GRCh38-aligned reads in locus windows +/- 1-2 Mb
plus all unmapped reads, so that differences reflect the recruitment rule, not the input):
  REF_ALN   reads whose PRIMARY alignment (winnowmap, GRCh38) overlaps the locus window W (= core +/- flank)
  REF_KMER  reads with >= F of their 31-mers in the GRCh38 W sequence (KMC 3 filter)
  POP_KMER  reads with >= F of their 31-mers in GRCh38 W + 40 HPRC release-2 haplotypes of W (sample's family excluded)
  ORACLE    reads with >= F of their 31-mers in the sample's own two truth haplotypes of W (headroom upper bound)
Each read set is assembled with hifiasm; hap1+hap2 contigs are aligned (minimap2 -cx asm5 --cs) to each truth
haplotype W; a truth haplotype is SPANNED if one alignment covers the whole projected core C (<=50 bp slack),
and RESOLVED if it is spanned at >= 99.9% identity across C.
Usage: 09_loci_pipeline.py --sample S --pool S.pool.bam --unmapped S.unmapped.bam --loci loci.tsv --proj truth_proj.tsv
       --truthdir DIR[,DIR] --popdir DIR --grch38dir DIR --gaps S_pat.gaps.bed,S_mat.gaps.bed --out OUTDIR [--threads 3] [--frac 0.5]
"""
import argparse, os, subprocess, tempfile, shutil, csv, re, time, json, glob
ap = argparse.ArgumentParser()
for k in ['sample', 'pool', 'unmapped', 'loci', 'proj', 'truthdir', 'popdir', 'grch38dir', 'gaps', 'out']: ap.add_argument('--' + k, required=True)
ap.add_argument('--threads', default='3'); ap.add_argument('--frac', default='0.5'); ap.add_argument('--only', default='')
ap.add_argument('--strategies', default='REF_ALN,REF_KMER,POP_KMER,ORACLE')
a = ap.parse_args()
T = a.threads
os.makedirs(a.out, exist_ok=True)
tmp = tempfile.mkdtemp(prefix='scratch.', dir=os.environ.get('TMPDIR', os.getcwd()))
def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, executable='/bin/bash', **kw).stdout
def fa_seqs(fn):
    d = {}; n = None
    for l in open(fn):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
cs_re = re.compile(r'(:\d+|\*[a-z][a-z]|[+\-][a-z]+)')
def aln_core(paf, c0, c1):
    """best single alignment over truth core [c0,c1): returns (covered_frac, identity_in_core, spanned, contig)"""
    best = (0.0, 0.0, False, '')
    for l in open(paf):
        f = l.rstrip('\n').split('\t')
        if len(f) < 12: continue
        ts, te = int(f[7]), int(f[8])
        ov = max(0, min(te, c1) - max(ts, c0))
        if ov <= 0: continue
        cs = next((x[5:] for x in f[12:] if x.startswith('cs:Z:')), '')
        pos = ts; m = e = 0
        for op in cs_re.findall(cs):
            c = op[0]
            if c == ':':
                n = int(op[1:]); a0, a1 = max(pos, c0), min(pos + n, c1)
                if a1 > a0: m += a1 - a0
                pos += n
            elif c == '*':
                if c0 <= pos < c1: e += 1
                pos += 1
            elif c == '+':
                if c0 <= pos < c1: e += len(op) - 1
            elif c == '-':
                n = len(op) - 1; a0, a1 = max(pos, c0), min(pos + n, c1)
                if a1 > a0: e += a1 - a0
                pos += n
        ident = m / max(1, m + e); frac = ov / (c1 - c0)
        spanned = ts <= c0 + 50 and te >= c1 - 50
        cand = (frac, ident, spanned, f[0])
        if (cand[2], cand[0], cand[1]) > (best[2], best[0], best[1]): best = cand
    return best
# ---- inputs
loci = [l for l in csv.DictReader(open(a.loci), delimiter='\t') if not a.only or l['locus'] in a.only.split(',')]
proj = {}
for l in csv.DictReader(open(a.proj), delimiter='\t'):
    if l['sample'] == a.sample: proj[(l['locus'], l['hap'], l['kind'])] = l
gaps = []
for g in a.gaps.split(','):
    for l in open(g):
        f = l.split(); gaps.append((f[0], int(f[1]), int(f[2])))
rows = []
t0 = time.time()
poolfq = os.path.join(tmp, 'pool.fq')
sh(f"(samtools fastq -F 0x900 {a.pool} 2>/dev/null; samtools fastq {a.unmapped} 2>/dev/null) | seqkit rmdup -n 2>/dev/null > {poolfq}")
pool_stats = sh(f"seqkit stats -T {poolfq}").split('\n')[1].split('\t')
pool_bp = int(pool_stats[4]); pool_n = int(pool_stats[3])
print(f'pool {pool_n} reads {pool_bp} bp in {time.time()-t0:.0f}s', flush=True)
try:
    for l in loci:
        loc = l['locus']; chrom = l['chrom']; ws, we = int(l['win_start']), int(l['win_end'])
        d = os.path.join(tmp, loc); os.makedirs(d, exist_ok=True)
        # truth haplotypes
        truth = {}
        for td in a.truthdir.split(','):
            for h in ('1', '2'):
                fn = os.path.join(td, loc, f'truth.{a.sample}.{h}.W.fa')
                if os.path.exists(fn): truth[h] = fn
        tw = os.path.join(d, 'truthW.fa'); open(tw, 'w').write(''.join(open(truth[h]).read() for h in sorted(truth)))
        gw = os.path.join(a.grch38dir, f'{loc}.W.fa')
        popw = os.path.join(a.popdir, loc, 'pop.W.fa')
        popref = os.path.join(d, 'popref.fa'); open(popref, 'w').write(open(gw).read() + open(popw).read())
        n_pop = sum(1 for x in open(popw) if x.startswith('>'))
        sets = {}
        # REF_ALN
        names = os.path.join(d, 'refaln.names')
        sh(f"samtools view -F 0x904 {a.pool} {chrom}:{ws+1}-{we} | cut -f1 | sort -u > {names}")
        sh(f"seqkit grep -f {names} {poolfq} > {d}/REF_ALN.fq")
        supp_only = int(sh(f"samtools view -f 0x800 {a.pool} {chrom}:{ws+1}-{we} | cut -f1 | sort -u | comm -23 - {names} | wc -l").strip())
        sets['REF_ALN'] = f'{d}/REF_ALN.fq'
        # k-mer recruitment
        kt = {}
        for strat, src in (('REF_KMER', gw), ('POP_KMER', popref), ('ORACLE', tw)):
            if strat not in a.strategies: continue
            db = os.path.join(d, strat); kd = os.path.join(d, strat + '_tmp'); os.makedirs(kd, exist_ok=True)
            sh(f"kmc -k31 -ci1 -cs65535 -fm -t{T} -m4 {src} {db} {kd} >/dev/null")
            s = time.time()
            sh(f"kmc_tools -t{T} filter {db} -ci1 {poolfq} -ci{a.frac} -cx1.0 {d}/{strat}.fq")
            kt[strat] = time.time() - s
            sets[strat] = f'{d}/{strat}.fq'
        names_of = {k: set(sh(f"awk 'NR%4==1' {v} | cut -c2- | cut -d' ' -f1").split()) for k, v in sets.items()}
        orc = names_of.get('ORACLE', set())
        for strat in [s for s in a.strategies.split(',') if s in sets]:
            fq = sets[strat]; st = sh(f"seqkit stats -T {fq}").split('\n')[1].split('\t')
            nreads, nbp = int(st[3]), int(st[4])
            row = dict(sample=a.sample, locus=loc, strategy=strat, n_reads=nreads, read_bp=nbp, pool_reads=pool_n,
                       n_pop_haps=n_pop, kmer_filter_s=round(kt.get(strat, 0), 2), supp_only_reads=supp_only if strat == 'REF_ALN' else '',
                       oracle_overlap=len(names_of[strat] & orc), not_in_oracle=len(names_of[strat] - orc),
                       oracle_missed=len(orc - names_of[strat]))
            # assemble
            ad = os.path.join(d, 'asm_' + strat); os.makedirs(ad, exist_ok=True)
            asmfa = os.path.join(ad, 'asm.fa'); open(asmfa, 'w').close(); row['asm_s'] = ''
            if nreads >= 5:
                try:
                    s = time.time()
                    sh(f"hifiasm -o {ad}/a -t{T} -f0 {fq} 2> {ad}/log")   # -f0: no 16-GB Bloom filter for tiny inputs
                    row['asm_s'] = round(time.time() - s, 1)
                    sh(f"cat {ad}/a.bp.hap1.p_ctg.gfa {ad}/a.bp.hap2.p_ctg.gfa | awk '/^S/{{print \">\"$2; print $3}}' > {asmfa}")
                except subprocess.CalledProcessError:
                    open(asmfa, 'w').close()
            seqs = fa_seqs(asmfa) if os.path.getsize(asmfa) else {}
            row['n_ctg'] = len(seqs); row['asm_bp'] = sum(len(s) for s in seqs.values())
            shutil.copy(asmfa, os.path.join(a.out, f'{a.sample}.{loc}.{strat}.asm.fa'))
            for h in sorted(truth):
                pw = proj.get((loc, h, 'W')); pc = proj.get((loc, h, 'C'))
                tseq = fa_seqs(truth[h]); tname = list(tseq)[0]; tlen = len(tseq[tname])
                if pc and pc.get('status') == 'ok' and pw:
                    if pw['strand'] == '+': c0, c1 = int(pc['tstart']) - int(pw['tstart']), int(pc['tend']) - int(pw['tstart'])
                    else: c0, c1 = int(pw['tend']) - int(pc['tend']), int(pw['tend']) - int(pc['tstart'])
                    csrc = 'projected'
                else:
                    fl = int(l['flank']); c0, c1 = fl, tlen - fl; csrc = 'W_minus_flank'
                c0, c1 = max(0, c0), min(tlen, c1)
                contig = pw['contig'] if pw else ''
                tg = sum(1 for g in gaps if g[0] == contig and g[1] < int(pw['tend']) and g[2] > int(pw['tstart'])) if pw else 0
                r = dict(row, hap=h, truth_W_len=tlen, core_len_truth=c1 - c0, core_src=csrc, truth_gaps_in_W=tg,
                         truth_W_status=pw['status'] if pw else 'missing')
                if seqs:
                    paf = os.path.join(ad, f'h{h}.paf')
                    sh(f"minimap2 -cx asm5 --cs -t{T} {truth[h]} {asmfa} > {paf} 2>/dev/null")
                    frac, ident, spanned, ctg = aln_core(paf, c0, c1)
                else: frac, ident, spanned, ctg = 0.0, 0.0, False, ''
                r.update(core_best_single_contig_frac=round(frac, 4), core_identity=round(ident, 6), spanned=int(spanned),
                         resolved=int(spanned and ident >= 0.999), best_contig=ctg)
                rows.append(r); print(json.dumps({k: r[k] for k in ('locus', 'strategy', 'hap', 'n_reads', 'n_ctg', 'spanned', 'resolved', 'core_identity', 'core_best_single_contig_frac')}), flush=True)
        shutil.rmtree(d, ignore_errors=True)
finally:
    keys = ['sample', 'locus', 'strategy', 'hap', 'n_reads', 'read_bp', 'pool_reads', 'n_pop_haps', 'kmer_filter_s', 'supp_only_reads',
            'oracle_overlap', 'not_in_oracle', 'oracle_missed', 'asm_s', 'n_ctg', 'asm_bp', 'truth_W_len', 'core_len_truth', 'core_src',
            'truth_gaps_in_W', 'truth_W_status', 'core_best_single_contig_frac', 'core_identity', 'spanned', 'resolved', 'best_contig']
    out = os.path.join(a.out, f'{a.sample}.loci.tsv')
    with open(out, 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows: f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
    with open(os.path.join(a.out, f'{a.sample}.pool.txt'), 'w') as f: f.write(f'pool_reads\t{pool_n}\npool_bp\t{pool_bp}\n')
    shutil.rmtree(tmp, ignore_errors=True)
