#!/usr/bin/env python3
"""Collate titration + loci results into L4_titration.tsv / L4_loci.tsv and draw two compressed figures.
Usage: 12_tables_figures.py <WORKDIR> <OUTDIR>"""
import sys, os, glob, json, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W, O = sys.argv[1], sys.argv[2]; os.makedirs(O, exist_ok=True)
raw = pd.read_csv(os.path.join(W, 'eval', 'titration_raw.tsv'), sep='\t')
raw[['sample', 'mode', 'cov']] = raw.label.str.split('|', expand=True)
rows = []
for (s, m, c), g in raw.groupby(['sample', 'mode', 'cov'], sort=False):
    h1 = g[g.file == 'hap1'].iloc[0]; h2 = g[g.file == 'hap2'].iloc[0]; b = g[g.file == 'hap1+hap2'].iloc[0]
    sub = os.path.join(W, 'sub', f'{s}.{m}.subsample.tsv'); ach = np.nan
    if os.path.exists(sub):
        t = pd.read_csv(sub, sep='\t', dtype={'target': str}); x = t[t.target == c.rstrip('x')]
        if len(x): ach = float(x.achieved_cov.iloc[0])
    tf = os.path.join(W, 'asm', f'{s}.{m}.{c}', 'time.txt'); wall = cpu = rss = np.nan
    if os.path.exists(tf):
        txt = open(tf).read()
        mm = re.search(r'Real time: ([\d.]+) sec; CPU: ([\d.]+) sec; Peak RSS: ([\d.]+) GB', txt)
        if mm: wall, cpu, rss = float(mm.group(1)), float(mm.group(2)), float(mm.group(3))
    rows.append(dict(sample=s, mode=m, target_cov=c, achieved_cov=ach,
                     n_ctg_hap1=h1.n_ctg, n_ctg_hap2=h2.n_ctg, NG50_hap1=h1.NG50, NG50_hap2=h2.NG50,
                     NG50_mean=(h1.NG50 + h2.NG50) / 2, NGA50_hap1=h1.NGA50, NGA50_hap2=h2.NGA50, NGA50_mean=(h1.NGA50 + h2.NGA50) / 2,
                     cov999_pat=b.cov999_pat, cov999_mat=b.cov999_mat, cov999_mean=b.cov999_mean, aln_pat=b.aln_pat, aln_mat=b.aln_mat,
                     QV_hap1=h1.QV, QV_hap2=h2.QV, switch_hap1=h1.switch_err, switch_hap2=h2.switch_err,
                     hamming_hap1=h1.hamming_err, hamming_hap2=h2.hamming_err,
                     hamming_truth_hap1=h1.get('hamming_truth'), hamming_truth_hap2=h2.get('hamming_truth'),
                     switch_truth_hap1=h1.get('switch_truth'), switch_truth_hap2=h2.get('switch_truth'),
                     wall_s=wall, cpu_s=cpu, max_rss_GB=rss))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(O, 'L4_titration.tsv'), sep='\t', index=False, float_format='%.6g')
# --- figure 1: titration
plt.rcParams.update({'font.size': 8})
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4))
sty = {('NA19240', 'hifi'): ('HiFi (NA19240)', '#1f77b4', 'o'), ('HG01258', 'hifi'): ('HiFi (HG01258)', '#17becf', 's'),
       ('HG01258', 'hifi_ul'): ('HiFi + ONT-R10 as --ul (HG01258)', '#2ca02c', '^'), ('HG01258', 'ont'): ('ONT-R10 only (HG01258)', '#d62728', 'D')}
for (s, m), (lab, col, mk) in sty.items():
    g = T[(T['sample'] == s) & (T['mode'] == m)].copy()
    if not len(g): continue
    g = g.sort_values('achieved_cov')
    if g.achieved_cov.isna().all(): continue
    x = g.achieved_cov
    axs[0, 0].plot(x, 100 * g.cov999_mean, marker=mk, color=col, label=lab, lw=1.2, ms=4)
    axs[0, 1].plot(x, g.NGA50_mean / 1e6, marker=mk, color=col, lw=1.2, ms=4)
    axs[1, 0].plot(x, (g.QV_hap1 + g.QV_hap2) / 2, marker=mk, color=col, lw=1.2, ms=4)
    axs[1, 1].plot(x, 100 * (g.switch_truth_hap1 + g.switch_truth_hap2) / 2, marker=mk, color=col, lw=1.2, ms=4)
axs[0, 0].set_ylabel('truth chr20 covered at >=99.9% id (%)'); axs[0, 1].set_ylabel('NGA50 (Mb, mean of hap1/hap2)')
axs[0, 1].set_yscale('log'); axs[1, 0].set_ylabel('QV (yak, Illumina chr20 k-mers)'); axs[1, 1].set_ylabel('phase switch rate, truth windows (%)')
tq = T[T['mode'] == 'truth']
if len(tq): axs[1, 0].axhline(float(((tq.QV_hap1 + tq.QV_hap2) / 2).mean()), ls='--', lw=0.8, color='grey', label='truth haplotypes (QV ceiling)'); axs[1, 0].legend(fontsize=6.5)
for ax in axs.flat: ax.set_xlabel('coverage of primary read type (x, chr20)'); ax.grid(alpha=.3)
axs[0, 0].legend(fontsize=6.5, loc='lower right')
fig.suptitle('L4: chr20 assembly vs coverage and read type (hifiasm 0.25.0; truth = HPRC r2 trio haplotypes)', fontsize=8.5)
fig.tight_layout(); fig.savefig(os.path.join(O, 'L4_fig1_titration.png'), dpi=170); plt.close(fig)
# --- loci
lf = sorted(glob.glob(os.path.join(W, 'loci_out', '*.loci.tsv')))
if lf:
    Lc = pd.concat([pd.read_csv(f, sep='\t') for f in lf])
    Lc.to_csv(os.path.join(O, 'L4_loci.tsv'), sep='\t', index=False)
    strats = ['REF_ALN', 'REF_KMER', 'POP_KMER', 'ORACLE']
    Lc['hapid'] = Lc['sample'] + '#' + Lc['hap'].astype(str)
    loci = list(dict.fromkeys(Lc.locus)); haps = sorted(Lc.hapid.unique())
    fig, axs = plt.subplots(1, len(strats), figsize=(7.4, 4.0), sharey=True)
    for k, s in enumerate(strats):
        M = np.full((len(loci), len(haps)), np.nan)
        for i, l in enumerate(loci):
            for j, h in enumerate(haps):
                x = Lc[(Lc.locus == l) & (Lc.hapid == h) & (Lc.strategy == s)]
                if len(x): M[i, j] = 2 * x.resolved.iloc[0] + (1 - x.resolved.iloc[0]) * x.spanned.iloc[0] if x.truth_gaps_in_W.iloc[0] == 0 else np.nan
        axs[k].imshow(M, cmap=matplotlib.colors.ListedColormap(['#f2f2f2', '#9ecae1', '#08519c']), vmin=0, vmax=2, aspect='auto')
        axs[k].set_xticks(range(len(haps))); axs[k].set_xticklabels([h.replace('#1', ' pat').replace('#2', ' mat') for h in haps], rotation=90, fontsize=6)
        axs[k].set_title(s.replace('_', '-'), fontsize=8)
        n_res = int(np.nansum(M == 2)); n_sp = int(np.nansum(M >= 1)); n = int(np.sum(~np.isnan(M)))
        axs[k].set_xlabel(f'resolved {n_res}/{n}\nspanned {n_sp}/{n}', fontsize=7)
    axs[0].set_yticks(range(len(loci))); axs[0].set_yticklabels(loci, fontsize=6.5)
    fig.suptitle('L4: hard loci, HiFi-only hifiasm per recruitment rule (dark: spanned at >=99.9%; light: spanned <99.9%; grey: not spanned)', fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(O, 'L4_fig2_loci.png'), dpi=170); plt.close(fig)
print(T.to_string())
