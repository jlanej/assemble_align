"""L1_figure.py -- L1_density_sensitivity.png from results/index_chm13.tsv, L1_results.tsv, results/curves.tsv.
usage: python L1_figure.py RESULTS_DIR L1_results.tsv OUT.png
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:   # house style helpers if available (figure-style skill); plain matplotlib otherwise
    from kernel import apply_figure_style, panel_letter   # noqa
except Exception:
    apply_figure_style = None
    def panel_letter(ax, letter, **kw):
        ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')

NAMES = {'mm2': 'minimap2 2.31 (random, hash64)', 'hashw': 'minimap2 + invertible width-adapted hash', 'rand_mix': 'random (strong hash)', 'miniception': 'miniception',
         'ddecyc': 'double decycling', 'greedy': 'GreedyMini (min over strands)', 'mod': 'mod-minimizer',
         'modoc': 'open-closed mod-minimizer', 'oc': 'open-closed minimizer', 'osync': 'open syncmer',
         'csync': 'closed syncmer'}
COL = {'mm2': '#000000', 'hashw': '#4d4d4d', 'rand_mix': '#7f7f7f', 'miniception': '#8c6d31', 'ddecyc': '#9467bd', 'greedy': '#bcbd22',
       'mod': '#1f77b4', 'modoc': '#d62728', 'oc': '#ff7f0e', 'osync': '#2ca02c', 'csync': '#17becf'}
MK = {'mm2': 'o', 'hashw': 'h', 'rand_mix': 'o', 'miniception': 's', 'ddecyc': 'D', 'greedy': 'P', 'mod': '^', 'modoc': 'v',
      'oc': '<', 'osync': 'X', 'csync': '*'}
PRESETS = ['sr', 'map-ont', 'map-hifi', 'asm20', 'lr:hqae']


def rtsv(p):
    L = open(p).read().strip().split('\n'); k = L[0].split('\t')
    return [dict(zip(k, l.split('\t'))) for l in L[1:]]


def main(res, restsv, out):
    if apply_figure_style:
        apply_figure_style(sizes=(8, 7, 6))
    idx = rtsv(os.path.join(res, 'index_chm13.tsv'))
    R = rtsv(restsv)
    cur = rtsv(os.path.join(res, 'curves.tsv')) if os.path.exists(os.path.join(res, 'curves.tsv')) else []
    fig, axs = plt.subplots(2, 3, figsize=(11.0, 6.6))
    # (a) density factor per preset
    ax = axs[0, 0]
    for j, p in enumerate(PRESETS):
        rows = [r for r in idx if r['preset'] == p]
        k, w = int(rows[0]['k']), int(rows[0]['w'])
        lb = np.ceil((w + k) / w) / (w + k) * (w + 1)
        ax.plot([j - 0.38, j + 0.38], [lb, lb], color='0.5', lw=1, ls='--')
        sch = [r['scheme'] for r in rows if r['scheme'] in NAMES]
        for i, r in enumerate([r for r in rows if r['scheme'] in NAMES]):
            x = j - 0.3 + 0.6 * i / max(1, len(sch) - 1)
            ax.scatter(x, float(r['dfactor']), color=COL[r['scheme']], marker=MK[r['scheme']], s=18, zorder=3,
                       edgecolor='k' if r['scheme'] in ('mm2', 'modoc') else 'none', linewidth=0.5)
    ax.axhline(2.0, color='0.75', lw=0.8, zorder=0)
    ax.set_xticks(range(len(PRESETS))); ax.set_xticklabels([f'{p}\n(k{rtsv_kw(idx, p)})' for p in PRESETS])
    ax.set_ylabel('density factor  d·(w+1)')
    ax.set_title('Lower-density schemes cut chr20 index entries', loc='left')
    ax.text(0.02, 0.02, 'grey line: random minimizer 2/(w+1)   dashed: forward-scheme lower bound',
            transform=ax.transAxes, fontsize=6, color='0.35', va='bottom')
    ax.set_ylim(0.95, 2.3)
    panel_letter(ax, 'a')
    # (b) entries change vs true-seed change (hap / asm reads)
    ax = axs[0, 1]
    for r in R:
        if r['scheme'] in ('mm2', 'mm2_1pick') or r['group'] not in ('hap', 'asm'):
            continue
        ax.scatter(float(r['entries_pct']), float(r['true_per_read_pct']), color=COL.get(r['scheme'], 'k'),
                   marker=MK.get(r['scheme'], 'o'), s=16, alpha=0.9)
    lim = [-32, 12]
    ax.plot(lim, lim, color='0.6', lw=0.8, ls=':')
    ax.axhline(0, color='0.8', lw=0.6); ax.axvline(0, color='0.8', lw=0.6)
    ax.set_xlabel('index entries vs minimap2 (%)'); ax.set_ylabel('true seeds per read vs minimap2 (%)')
    ax.set_title('True seeds track index size (dotted: y = x)', loc='left')
    panel_letter(ax, 'b')
    # (c) longest seedless gap p95 change vs entries change
    ax = axs[0, 2]
    for r in R:
        if r['scheme'] in ('mm2', 'mm2_1pick') or r['group'] not in ('hap:uniq', 'asm:uniq'):
            continue
        ax.scatter(float(r['entries_pct']), float(r['gap_p95_pct']), color=COL.get(r['scheme'], 'k'),
                   marker=MK.get(r['scheme'], 'o'), s=16)
    ax.axhline(0, color='0.8', lw=0.6); ax.axvline(0, color='0.8', lw=0.6)
    ax.set_xlabel('index entries vs minimap2 (%)'); ax.set_ylabel('longest-gap p95 vs minimap2 (%)')
    ax.set_title('Longest seedless gap, unique-sequence reads', loc='left')
    panel_letter(ax, 'c')
    # (d-f) equal-density curves
    for ax, p, lab, letter in [(axs[1, 0], 'sr', 'Illumina 2x150, NA19240 reads, k=21', 'd'),
                               (axs[1, 1], 'map-ont', 'ONT Q20 model, NA19240 reads, k=15', 'e'),
                               (axs[1, 2], 'map-hifi', 'HiFi Q30 model, NA19240 reads, k=19', 'f')]:
        for sch in ['mm2', 'mod', 'modoc', 'oc', 'miniception', 'osync', 'csync']:
            pts = sorted([(float(c['density']), float(c['true_per_kb'])) for c in cur if c['preset'] == p and c['scheme'] == sch])
            if not pts:
                continue
            x, y = zip(*pts)
            ax.plot(x, y, marker=MK[sch], ms=3.5, lw=1.2 if sch in ('mm2', 'modoc') else 0.8, color=COL[sch],
                    label=NAMES[sch])
        k, w = E_KW[p]
        ax.axvline([float(r['density']) for r in idx if r['preset'] == p and r['scheme'] == 'mm2'][0], color='0.8', lw=0.8, ls='--')
        ax.set_xlabel('density (fraction of k-mers indexed)'); ax.set_ylabel('true-origin seeds per kb of read')
        ax.set_title(lab, loc='left')
        if letter == 'd':
            ax.text(0.172, 0.04, 'minimap2 preset density', transform=ax.get_xaxis_transform(), fontsize=6, color='0.4')
        panel_letter(ax, letter)
    h, l = axs[1, 0].get_legend_handles_labels()
    extra = [plt.Line2D([], [], color=COL[s], marker=MK[s], ls='none', label=NAMES[s]) for s in ('rand_mix', 'hashw', 'ddecyc', 'greedy')]
    fig.legend(h + extra, l + [e.get_label() for e in extra], loc='lower center', ncol=5, frameon=False, fontsize=6.5,
               bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, dpi=170)
    print('saved', out)


E_KW = {'sr': (21, 11), 'map-ont': (15, 10), 'map-hifi': (19, 19), 'asm20': (19, 10), 'lr:hqae': (25, 51)}


def rtsv_kw(idx, p):
    k, w = E_KW[p]
    return f'{k} w{w}'


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
