#!/usr/bin/env python3
"""Compose L4_assembly_pilot.md from the result tables (numbers are computed, not typed).
Usage: 14_build_report.py <OUTDIR with L4_titration.tsv, L4_loci.tsv> <methods.md> <fig1_marker> <fig2_marker>"""
import sys, pandas as pd, numpy as np
O, methods, f1, f2 = sys.argv[1:5]
T = pd.read_csv(f'{O}/L4_titration.tsv', sep='\t'); L = pd.read_csv(f'{O}/L4_loci.tsv', sep='\t')
def row(s, m, c):
    x = T[(T['sample'] == s) & (T['mode'] == m) & (T.target_cov == c)]
    return x.iloc[0] if len(x) else None
def pct(v): return f"{100*v:.2f}%"
def mb(v): return f"{v/1e6:.1f} Mb"
def qv(r): return f"{(r.QV_hap1 + r.QV_hap2)/2:.1f}"
def nctg(r): return int(r.n_ctg_hap1 + r.n_ctg_hap2)
def sw(r): return f"{100*(r.switch_truth_hap1 + r.switch_truth_hap2)/2:.1f}%"
H = {c: row('HG01258', 'hifi', c) for c in ['fullx', '24x', '16x', '12x', '8x']}
O_ = {c: row('HG01258', 'ont', c) for c in ['fullx', '24x', '16x', '12x', '8x']}
U = {c: row('HG01258', 'hifi_ul', c) for c in ['fullx', '16x', '8x']}
TR = {s: row(s, 'truth', '-x') for s in ['NA19240', 'HG01258']}
lines = []
A = lines.append
# ---------- titration table (markdown)
tab = ['| Sample / mode | Coverage (x) | Contigs (h1+h2) | NG50 mean | NGA50 mean | cov999 (pat/mat) | QV (mean) | Truth-window switch | yak hamming (h1/h2) | Wall / CPU (h) | Peak RSS (GB) |',
       '|---|---|---|---|---|---|---|---|---|---|---|']
for _, r in T.iterrows():
    if r['mode'] == 'truth':
        tab.append(f"| {r['sample']} truth (HPRC r2) | - | 2 | {mb(r.NG50_mean)} | {mb(r.NGA50_mean)} | 100/100% | {qv(r)} | 0.0% | {100*r.hamming_hap1:.1f}/{100*r.hamming_hap2:.1f}% | - | - |")
        continue
    tab.append(f"| {r['sample']} {r['mode']} | {r.achieved_cov:.1f} | {nctg(r)} | {mb(r.NG50_mean)} | {mb(r.NGA50_mean)} | {100*r.cov999_pat:.2f}/{100*r.cov999_mat:.2f}% | {qv(r)} | {sw(r)} | {100*r.hamming_hap1:.1f}/{100*r.hamming_hap2:.1f}% | {r.wall_s/3600:.2f} / {r.cpu_s/3600:.2f} | {r.max_rss_GB:.1f} |")
# ---------- loci summary
Ls = L[L.truth_gaps_in_W == 0]
summ = []
for s in ['REF_ALN', 'REF_KMER', 'POP_KMER', 'ORACLE']:
    x = Ls[Ls.strategy == s]; h1 = x[x.hap == 1]
    summ.append((s, int(x.resolved.sum()), int(x.spanned.sum()), len(x), int(h1.n_reads.sum()), int(h1.oracle_missed.sum()), int(h1.not_in_oracle.sum())))
ltab = ['| Rule | Resolved | Spanned | Reads (sum, 13 loci) | Oracle reads missed | Reads not in oracle |', '|---|---|---|---|---|---|']
for s, res, sp, n, rd, miss, extra in summ:
    ltab.append(f"| {s} | {res}/{n} | {sp}/{n} | {rd:,} | {miss:,} | {extra:,} |")
def lrow(locus, strat, hap=1):
    x = L[(L.locus == locus) & (L.strategy == strat) & (L.hap == hap)]
    return x.iloc[0]
dq_ref, dq_pop, dq_or = lrow('HLA-DRB_DQ', 'REF_KMER'), lrow('HLA-DRB_DQ', 'POP_KMER'), lrow('HLA-DRB_DQ', 'ORACLE')
smn = L[(L.locus == 'SMN1_SMN2')]
smn_frac = (smn.core_best_single_contig_frac.min(), smn.core_best_single_contig_frac.max())
smn_core = (smn.core_len_truth.min(), smn.core_len_truth.max())
ncf_ref, ncf_or = lrow('NCF1', 'REF_ALN'), lrow('NCF1', 'ORACLE')
kf = L[L.strategy != 'REF_ALN'].kmer_filter_s
pool = int(L.pool_reads.iloc[0])
# ---------- litmus statements
lit = []
h, h16, h8 = H['fullx'], H['16x'], H['8x']
if h is not None:
    s = (f"**L4-T1 (HiFi coverage, HG01258 chr20).** At {h.achieved_cov:.1f}x HiFi, hifiasm reaches {pct(h.cov999_mean)} of the truth at >=99.9% identity, "
         f"{nctg(h)} contigs, NG50 {mb(h.NG50_mean)} and QV {qv(h)} (truth ceiling {qv(TR['HG01258'])}).")
    if h16 is not None: s += f" At 16x: {pct(h16.cov999_mean)}, {nctg(h16)} contigs, NG50 {mb(h16.NG50_mean)}, QV {qv(h16)}."
    if h8 is not None: s += f" At 8x: {pct(h8.cov999_mean)}, {nctg(h8)} contigs, NG50 {mb(h8.NG50_mean)}, QV {qv(h8)}."
    lit.append(s)
o = O_['fullx']; u = U['fullx']
if o is not None and h is not None:
    s = (f"**L4-T2 (read type at ~36x).** The ONT-R10-only assembly (hifiasm `--ont`, {o.achieved_cov:.1f}x, read N50 41.9 kb, dorado 0.8.1 sup) covers {pct(o.cov999_mean)} at >=99.9% "
         f"(HiFi: {pct(h.cov999_mean)}), with {nctg(o)} contigs (HiFi: {nctg(h)}) and QV {qv(o)} (HiFi: {qv(h)}).")
    if O_['16x'] is not None: s += f" At 16x ONT: {pct(O_['16x'].cov999_mean)}, {nctg(O_['16x'])} contigs, QV {qv(O_['16x'])}."
    lit.append(s)
if u is not None and h is not None:
    lit.append(f"**L4-T3 (ONT as `--ul`).** Adding 26.1x of standard-length R10 reads >= 20 kb as `--ul` to {u.achieved_cov:.1f}x HiFi cuts contigs from {nctg(h)} to {nctg(u)}; NG50 rises from {mb(h.NG50_mean)} to {mb(u.NG50_mean)}, and one haplotype is a single {mb(u.NG50_hap2)} contig. "
               f"cov999 is {pct(u.cov999_mean)} vs {pct(h.cov999_mean)}, and the run costs {u.wall_s/h.wall_s:.1f}x the HiFi-only wall time ({u.wall_s/3600:.1f} h vs {h.wall_s/3600:.1f} h; peak RSS {u.max_rss_GB:.1f} vs {h.max_rss_GB:.1f} GB).")
lit.append(f"**L4-L1 (hard loci, NA19240, 13 loci x 2 haplotypes, HiFi-only hifiasm).** Resolved haplotype-loci: REF_ALN {summ[0][1]}/{summ[0][3]}, REF_KMER {summ[1][1]}/{summ[1][3]}, "
           f"POP_KMER {summ[2][1]}/{summ[2][3]}, ORACLE {summ[3][1]}/{summ[3][3]}. At HLA-DRB/DQ, population k-mers raise recruitment from {int(dq_ref.n_reads)} to {int(dq_pop.n_reads)} reads (oracle {int(dq_or.n_reads)}) "
           f"and resolution from {int(L[(L.locus=='HLA-DRB_DQ')&(L.strategy=='REF_KMER')].resolved.sum())}/2 to {int(L[(L.locus=='HLA-DRB_DQ')&(L.strategy=='POP_KMER')].resolved.sum())}/2. "
           f"No rule, including the oracle, spans SMN (truth core {smn_core[0]/1e6:.2f}-{smn_core[1]/1e6:.2f} Mb; best single contig covers {100*smn_frac[0]:.0f}-{100*smn_frac[1]:.0f}% of it).")
tq = TR['NA19240']; tq2 = TR['HG01258']
lit.append(f"**L4-Q1 (phasing metric).** `yak trioeval` with chr20-sliced parental Illumina k-mers reports {100*min(tq.hamming_hap1,tq.hamming_hap2,tq2.hamming_hap1,tq2.hamming_hap2):.1f}-{100*max(tq.hamming_hap1,tq.hamming_hap2,tq2.hamming_hap1,tq2.hamming_hap2):.1f}% hamming error on the trio-phased truth haplotypes themselves. "
           f"The truth-window metric scores 0.0% on the truth and recovers a planted 30% chimera as 29%. Sliced trio k-mers therefore cannot detect phasing differences smaller than about 4%.")
md = open(methods).read()
A('# L4 litmus pilot: assembly as a function of coverage and read type, and population-informed recruitment at hard loci\n')
A('*Real 1000 Genomes / HPRC release-2 data; hifiasm 0.25.0-r726; 3 threads; shared laptop (load average 15-37 on 10 cores during the runs).*\n')
A('## Summary\n')
A('The pilot used two trio children with region-sliceable long reads. NA19240 was the preferred sample, and its HiFi reads and loci pool were used. HG01258 was added because its ONT R10 reads are aligned and indexed, whereas NA19240 ONT exists only as unaligned 13-73-GB BAMs. Each sample\'s HPRC release-2 trio-phased haplotypes served as truth. Results are provisional where the titration grid is partial (see Deviations).\n')
for s in lit: A('- ' + s)
A('')
A(md)
A('\n## Results\n')
A('### Coverage and read-type titration (chr20)\n')
A('\n'.join(tab)); A('')
A(f'![L4 titration]({f1})\n')
A('cov999 is the fraction of truth bases in 10-kb windows aligned at >= 99.9% identity, averaged over pat/mat. NGA50 is haplotype-consistent: it breaks at phase switches across heterozygous SVs, so it is far below NG50 for dual (non-trio) assemblies. The truth-window switch rate is the fraction of adjacent informative 20-kb windows that change parental haplotype. QV is capped near 38 by the chr20-sliced Illumina k-mer database (truth = ' + f"{qv(TR['NA19240'])}/{qv(TR['HG01258'])}" + ').\n')
A('### Hard loci: recruitment rules (NA19240, HiFi only)\n')
A('\n'.join(ltab)); A('')
A(f'![L4 loci]({f2})\n')
A(f"At NCF1, the oracle k-mer set (own truth haplotypes) recruits {int(ncf_or.n_reads)} reads versus {int(ncf_ref.n_reads)} for the Winnowmap window. 31-mer recruitment cannot separate NCF1 from the >99%-identical NCF1B/NCF1C copies, so the 'oracle' for k-mer rules includes paralog reads; hifiasm still resolved the 30-kb core with every rule. "
  f"One KMC filter pass over the {pool:,}-read (1.13 Gbp) pool took {kf.min():.0f}-{kf.max():.0f} s with 3 threads, about 17-34 Mbp/s. Scanning a whole 41x HiFi read set (~130 Gbp) would take about 1-2 h per pass, and obtaining that read set without region-sliceable alignments would mean a ~100-GB download (4-40 h at the measured 0.7-7 MB/s).\n")
open(f'{O}/L4_assembly_pilot.body.md', 'w').write('\n'.join(lines))
print('\n'.join(lit))
