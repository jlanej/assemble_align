#!/usr/bin/env python3
"""l5_figure.py - one summary figure for the L5 MAPQ-calibration pilot.
(a,b) empirical error vs reported MAPQ inside / outside SDs (GRCh38, CHM13);
(c) ECE inside SDs for raw MAPQ, MAPQ-only isotonic lookup, and the GBM + isotonic recalibrator;
(d) fraction of correctly placed in-SD reads retained at <=1e-3 error: raw MAPQ vs recalibrated.
usage: l5_figure.py OUTDIR (containing L5_calibration.tsv, L5_recalibrator_metrics.tsv) FIG.png"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUN_ORDER = ['bwa_sr', 'mm2_sr', 'mm2_hifi', 'mm2_lrhq', 'mm2_ont']
RUN_LABEL = {'bwa_sr': 'Illumina, bwa mem', 'mm2_sr': 'Illumina, minimap2 sr', 'mm2_hifi': 'HiFi, minimap2 map-hifi',
             'mm2_lrhq': 'ONT, minimap2 lr:hq', 'mm2_ont': 'ONT, minimap2 map-ont'}
COL = {'bwa_sr': '#1b6ca8', 'mm2_sr': '#6baed6', 'mm2_hifi': '#7b3294', 'mm2_lrhq': '#d95f02', 'mm2_ont': '#fdb863'}
BIN_X = {'0': 0, '1-3': 2, '4-9': 6.5, '10-19': 14.5, '20-29': 24.5, '30-39': 34.5, '40-49': 44.5, '50-59': 54.5, '60': 60}


def split_run(r):
    ref, run = r.split('.', 1)
    return ref, run


def plot(outdir, fig_path, dpi=170):
    cal = pd.read_csv(f'{outdir}/L5_calibration.tsv', sep='\t')
    if 'binned_by' in cal.columns:
        cal = cal[cal.binned_by == 'reported_mapq']
    met = pd.read_csv(f'{outdir}/L5_recalibrator_metrics.tsv', sep='\t')
    cal[['ref', 'cfg']] = cal.run.apply(lambda r: pd.Series(split_run(r)))
    met[['ref', 'cfg']] = met.run.apply(lambda r: pd.Series(split_run(r)))
    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(1, 2, wspace=0.08, left=0.09, right=0.98, top=0.93, bottom=0.585)
    gs2 = fig.add_gridspec(1, 2, wspace=0.08, left=0.22, right=0.98, top=0.435, bottom=0.08)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    xs = np.linspace(0, 60, 61)
    ttls = {}
    for ref in ['grch38', 'chm13']:
        v = cal[(cal.ref == ref) & (cal.stratum == 'in_SD') & (cal.mapq_bin.astype(str) == '60')].w_error
        ttls[ref] = '%s: in-SD MAPQ 60 reads err %.1f-%.1f%%' % ('GRCh38' if ref == 'grch38' else 'CHM13', 100 * v.min(), 100 * v.max())
    for ax, ref in zip(axes, ['grch38', 'chm13']):
        ax.plot(xs, 10 ** (-xs / 10), color='0.55', lw=0.8, ls=':', zorder=1)
        for cfg in RUN_ORDER:
            for stratum, ls, mk in (('in_SD', '-', 'o'), ('outside_SD', '--', 's')):
                c = cal[(cal.ref == ref) & (cal.cfg == cfg) & (cal.stratum == stratum)].copy()
                if c.empty:
                    continue
                c['x'] = c.mapq_bin.astype(str).map(BIN_X)
                c = c.sort_values('x')
                y = c.w_error.values.copy()
                zero = y <= 0
                y[zero] = c.err_ci_hi.values[zero]            # zero observed errors -> plot 95% upper bound
                ax.plot(c.x, y, ls=ls, color=COL[cfg], lw=1.1, zorder=2)
                ax.scatter(c.x[~zero], y[~zero], marker=mk, s=12, color=COL[cfg], zorder=3,
                           edgecolor='white', linewidth=0.3)
                ax.scatter(c.x[zero], y[zero], marker='v', s=14, facecolor='none', edgecolor=COL[cfg], zorder=3,
                           linewidth=0.8)
        ax.set_yscale('log'); ax.set_ylim(3e-7, 1.5); ax.set_xlim(-2, 63)
        ax.set_xlabel('Reported MAPQ (binned)')
        ax.set_title(ttls[ref], loc='left', fontsize=8)
    axes[0].set_ylabel('Empirical error rate')
    axes[1].tick_params(labelleft=False)
    axes[0].text(1, 4e-7, 'dotted: perfect calibration', fontsize=6, color='0.4')
    # legend (configs + strata)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=COL[c], lw=1.4, label=RUN_LABEL[c]) for c in RUN_ORDER]
    h += [Line2D([], [], color='0.3', ls='-', marker='o', ms=3, label='truth >=50% in SD'),
          Line2D([], [], color='0.3', ls='--', marker='s', ms=3, label='outside SD'),
          Line2D([], [], color='0.3', ls='none', marker='v', mfc='none', ms=4, label='0 errors: 95% upper bound')]
    axes[1].legend(handles=h, fontsize=5.5, frameon=False, loc='lower left', bbox_to_anchor=(0.0, 0.04), ncol=1)
    # (c) ECE inside SD
    ax = fig.add_subplot(gs2[0, 0])
    m = met[met.stratum == 'in_SD']
    rows = [(ref, cfg) for ref in ['grch38', 'chm13'] for cfg in RUN_ORDER
            if not m[(m.ref == ref) & (m.cfg == cfg)].empty]
    ylab = []
    for i, (ref, cfg) in enumerate(rows):
        mm = m[(m.ref == ref) & (m.cfg == cfg)].set_index('predictor')
        vals = [mm.loc[p, 'ECE'] for p in ('raw', 'iso', 'gbmcal' if 'gbmcal' in mm.index else 'gbm')]
        ax.plot(vals, [i] * 3, color='0.75', lw=0.8, zorder=1)
        ax.scatter(vals[0], i, marker='o', s=18, color='0.35', zorder=2)
        ax.scatter(vals[1], i, marker='D', s=14, facecolor='white', edgecolor='0.35', zorder=2)
        ax.scatter(vals[2], i, marker='o', s=22, color=COL[cfg], zorder=3)
        ylab.append(('GRCh38 ' if ref == 'grch38' else 'CHM13 ') + RUN_LABEL[cfg].split(', ')[1])
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(ylab, fontsize=5.5); ax.invert_yaxis()
    ax.set_xscale('log'); ax.set_xlabel('In-SD expected calibration error')
    ax.text(0.98, 0.02, 'lower = better', transform=ax.transAxes, ha='right', va='bottom', fontsize=6, color='0.4')
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
    tk = [0.01, 0.02, 0.05, 0.1, 0.2]
    ax.xaxis.set_major_locator(FixedLocator(tk)); ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in tk]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title('Recalibration shrinks in-SD calibration error', loc='left', fontsize=8)
    shared = [Line2D([], [], color='0.35', marker='o', ls='none', ms=4, label='reported MAPQ'),
              Line2D([], [], color='0.35', marker='D', mfc='white', ls='none', ms=3.5, label='MAPQ-only lookup (c)'),
              Line2D([], [], color='0.1', marker='o', ls='none', ms=4.5, label='recalibrator\n(colour = aligner;\nc: GBM + isotonic,\nd: GBM score)')]
    # (d) retained at 1e-3, in SD
    ax = fig.add_subplot(gs2[0, 1])
    col = 'frac_correct_retained_err1e-03'
    for i, (ref, cfg) in enumerate(rows):
        mm = m[(m.ref == ref) & (m.cfg == cfg)].set_index('predictor')
        a_, b_ = mm.loc['raw', col], mm.loc['gbm', col]
        ax.plot([a_, b_], [i, i], color='0.75', lw=0.8, zorder=1)
        ax.scatter(a_, i, marker='o', s=18, color='0.35', zorder=2)
        ax.scatter(b_, i, marker='o', s=22, color=COL[cfg], zorder=3)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([]); ax.invert_yaxis()
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel('Correct in-SD reads kept at \u22640.1% error')
    ax.legend(handles=shared, fontsize=5.5, frameon=False, loc='upper left', bbox_to_anchor=(0.64, 0.99),
              handletextpad=0.3, borderaxespad=0.0)
    ax.text(0.98, 0.02, 'higher = better', transform=ax.transAxes, ha='right', va='bottom', fontsize=6, color='0.4')
    ax.set_title('No MAPQ cut-off reaches 0.1% error in SDs', loc='left', fontsize=8)
    for a, l in zip(fig.axes, 'abcd'):
        a.text(-0.02, 1.07, l, transform=a.transAxes, fontsize=10, fontweight='bold', ha='right', va='bottom')
    fig.savefig(fig_path, dpi=dpi)
    return fig


if __name__ == '__main__':
    plot(sys.argv[1], sys.argv[2])
