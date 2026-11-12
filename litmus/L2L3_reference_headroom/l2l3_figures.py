#!/usr/bin/env python3
"""L2/L3 pilot figures (compressed PNG, 150 dpi).

  l2l3_figures.py --results results/ --outdir figures/
Fig 1  L2L3_fig1_reference_headroom.png : per technology (columns), per reference (bars):
       fraction of reads with MAPQ<20, small-variant edit rate, and GRCh38-REF allele
       fraction at truth heterozygous SNVs (MAPQ>=20) by stratum.
Fig 2  L2L3_fig2_longread_residuals.png : residual HiFi/ONT errors versus the sample's own
       assembly by sequence context, and reported-vs-empirical base-quality calibration.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REFS = ['GRCh38', 'CHM13', 'pers2p', 'pers1KG', 'oracle']
REF_LABEL = {'GRCh38': 'GRCh38', 'CHM13': 'CHM13', 'pers2p': 'GRCh38 + own calls\n(unphased)',
             'pers1KG': 'GRCh38 + 1kGP\nphased calls', 'oracle': 'own assembly\n(oracle)'}
REF_SHORT = {'GRCh38': 'GRCh38', 'CHM13': 'CHM13', 'pers2p': '+own calls', 'pers1KG': '+1kGP phased',
             'oracle': 'own asm'}
COL = {'GRCh38': '#7f7f7f', 'CHM13': '#bdbdbd', 'pers2p': '#9ecae1', 'pers1KG': '#3182bd', 'oracle': '#08306b'}
TECH = {'ilmn': 'Illumina 2x150', 'hifi': 'PacBio HiFi (DeepConsensus)', 'ont': 'ONT R10.4.1 (UL, sup)'}
STRATA = ['ALL', 'SD', 'TR', 'other']


def style():
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6,
                         'ytick.labelsize': 6, 'legend.fontsize': 7, 'axes.spines.top': False,
                         'axes.spines.right': False, 'xtick.direction': 'out', 'ytick.direction': 'out',
                         'legend.frameon': False, 'axes.titlelocation': 'left', 'savefig.dpi': 150})


def fig1(res, out):
    m = pd.read_csv(os.path.join(res, 'L2L3_metrics.tsv'), sep='\t')
    ab = pd.read_csv(os.path.join(res, 'L2L3_allele_balance.tsv'), sep='\t')
    ab = ab[ab.mapq_min.astype(str) == '20']
    if 'siteset' in ab.columns:
        ab = ab[ab.siteset == 'clean']
    techs = [t for t in ['ilmn', 'hifi', 'ont'] if t in set(m.tech)]
    fig, axes = plt.subplots(3, len(techs), figsize=(7.2, 7.6), squeeze=False)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.95, bottom=0.115, hspace=0.62, wspace=0.3)
    nread = {t: int(m[(m.tech == t) & (m.ref == 'GRCh38')].n_reads.iloc[0]) for t in techs}
    fmt = lambda n: ('%.1fM' % (n / 1e6)) if n >= 1e6 else ('%.1fk' % (n / 1e3))
    handles = None
    for j, t in enumerate(techs):
        mt = m[m.tech == t].set_index('ref')
        refs = [r for r in REFS if r in mt.index]
        x = np.arange(len(refs))
        for i, (col, ylab, f) in enumerate([('frac_mapq_lt20', 'reads with MAPQ<20\nor unmapped (%)', '%.2f'),
                                            ('edit_rate_small_pooled', 'edit rate, indels\n<50 bp (%)', '%.3f')]):
            ax = axes[i, j]
            v = [100 * mt.loc[r, col] for r in refs]
            ax.bar(x, v, color=[COL[r] for r in refs], width=0.7)
            for xi, vi in zip(x, v):
                ax.text(xi, vi, f % vi, ha='center', va='bottom', fontsize=6)
            ax.set_xticks(x, [REF_SHORT[r] for r in refs], rotation=30, ha='right')
            ax.margins(y=0.16)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title('%s\n(n = %s reads)' % (TECH[t], fmt(nread[t])), loc='left')
        ax = axes[2, j]
        at = ab[ab.tech == t]
        w = 0.8 / len(refs)
        bars = []
        for k, r in enumerate(refs):
            y = [at[(at.ref == r) & (at.stratum == s)].ref_frac_pooled.values for s in STRATA]
            y = [yy[0] if len(yy) else np.nan for yy in y]
            bars.append(ax.bar(np.arange(len(STRATA)) + (k - (len(refs) - 1) / 2) * w, y, width=w, color=COL[r],
                               label=REF_LABEL[r].replace('\n', ' ')))
        ax.axhline(0.5, color='k', lw=0.6, ls='--')
        ax.set_xticks(np.arange(len(STRATA)), ['all (no CEN)', 'seg. dup.', 'tandem rep.', 'other'], rotation=30, ha='right')
        vals = at[at.ref.isin(refs) & at.stratum.isin(STRATA)].ref_frac_pooled.values
        ax.set_ylim(min(0.47, np.nanmin(vals) - 0.01), max(0.53, np.nanmax(vals) + 0.01))
        if j == 0:
            ax.set_ylabel('GRCh38-REF allele fraction\n(clean het SNVs, MAPQ>=20)')
        handles = bars
    fig.legend(handles, [h.get_label() for h in handles], loc='lower center', ncol=5, fontsize=7,
               handlelength=1.0, columnspacing=1.0, bbox_to_anchor=(0.53, 0.0))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig2(res, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    order = ['HP1', 'HP2-3', 'STR', 'HP4-5', 'HP6-7', 'HP8-10', 'HP11-15', 'HP16+']
    cols = {'hifi': '#d95f02', 'ont': '#1b9e77'}
    lab = {'hifi': 'HiFi (DeepConsensus)', 'ont': 'ONT R10.4.1 sup'}
    ax = axes[0]
    for t, off in (('hifi', -0.2), ('ont', 0.2)):
        f = os.path.join(res, 'L2L3_lr_error_context_%s.tsv' % t)
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep='\t').set_index('context')
        o = [c for c in order if c in d.index]
        ax.scatter(np.arange(len(o)) + off, d.loc[o, 'err_per_kb'], color=cols[t], s=14, label=lab[t], zorder=3)
        ax.plot(np.arange(len(o)) + off, d.loc[o, 'err_per_kb'], color=cols[t], lw=0.6)
    ax.set_yscale('log')
    ax.set_xticks(np.arange(len(order)), ['non-HP', 'HP 2-3', 'STR', 'HP 4-5', 'HP 6-7', 'HP 8-10', 'HP 11-15', 'HP 16+'], rotation=35, ha='right')
    ax.set_ylabel('read errors per kb of assembly\n(mismatch + indel bases)')
    ax.set_title('Errors vs own assembly rise with HP length')
    ax.legend(loc='upper left')
    ax = axes[1]
    for t in ('hifi', 'ont'):
        f = os.path.join(res, 'L2L3_lr_qv_calibration_%s.tsv' % t)
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep='\t')
        d = d[d.bases >= 10000]
        ax.plot(d.reported_qv, d.empirical_qv, marker='o', ms=2.5, lw=0.8, color=cols[t], label=lab[t])
    lim = [0, 60]
    ax.plot(lim, lim, color='k', lw=0.6, ls='--')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('reported base QV')
    ax.set_ylabel('empirical QV vs own assembly')
    ax.set_title('Reported vs empirical base QV')
    ax.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    fig1(a.results, os.path.join(a.outdir, 'L2L3_fig1_reference_headroom.png'))
    fig2(a.results, os.path.join(a.outdir, 'L2L3_fig2_longread_residuals.png'))


if __name__ == '__main__':
    main()
