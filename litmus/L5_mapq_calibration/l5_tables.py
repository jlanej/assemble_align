#!/usr/bin/env python3
"""l5_tables.py - render the L5 markdown result tables from the TSV outputs (no hand-typed numbers).
usage: l5_tables.py RECAL_DIR RECAL_STRICT_DIR > tables.md"""
import sys
import numpy as np
import pandas as pd

RUNS = ['grch38.bwa_sr', 'grch38.mm2_sr', 'grch38.mm2_hifi', 'grch38.mm2_lrhq', 'grch38.mm2_ont',
        'chm13.bwa_sr', 'chm13.mm2_sr', 'chm13.mm2_hifi', 'chm13.mm2_lrhq', 'chm13.mm2_ont']
LAB = {'bwa_sr': 'Illumina · bwa mem', 'mm2_sr': 'Illumina · minimap2 sr', 'mm2_hifi': 'HiFi · minimap2 map-hifi',
       'mm2_lrhq': 'ONT · minimap2 lr:hq', 'mm2_ont': 'ONT · minimap2 map-ont'}


def lab(run):
    ref, cfg = run.split('.', 1)
    return ('GRCh38' if ref == 'grch38' else 'CHM13') + ' · ' + LAB[cfg]


def pct(x, d=2):
    return '—' if pd.isna(x) else ('%.' + str(d) + 'f%%') % (100 * x)


def sci(x):
    if pd.isna(x):
        return '—'
    if x == 0:
        return '0'
    return '%.2g' % x


def phred(e):
    return '—' if (pd.isna(e) or e <= 0) else '%.0f' % (-10 * np.log10(e))


def md(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    return '\n'.join(out)


def main(rd, rs):
    m = pd.read_csv(f'{rd}/L5_recalibrator_metrics.tsv', sep='\t')
    ms = pd.read_csv(f'{rs}/L5_recalibrator_metrics.tsv', sep='\t')
    dg = pd.read_csv(f'{rd}/L5_q60_error_diagnosis_summary.tsv', sep='\t')
    g = lambda df, run, st, p, c: df[(df.run == run) & (df.stratum == st) & (df.predictor == p)][c].iloc[0]
    # Table 1: raw calibration at MAPQ 0 / 60, in vs outside SD
    rows = []
    for run in RUNS:
        rows.append({'Run': lab(run),
                     'Evaluable reads': f"{g(m, run, 'all', 'raw', 'n_reads'):,}",
                     'MAPQ 0: n (error)': f"{g(m, run, 'all', 'raw', 'n_mapq0'):,} ({pct(g(m, run, 'all', 'raw', 'w_err_mapq0'), 1)})",
                     'MAPQ 60: n, wrong': f"{g(m, run, 'all', 'raw', 'n_mapq60'):,}, {g(m, run, 'all', 'raw', 'n_wrong_mapq60'):,}",
                     'MAPQ 60 error, outside SD': sci(g(m, run, 'outside_SD', 'raw', 'w_err_mapq60')),
                     'MAPQ 60 in SD: wrong / n': f"{g(m, run, 'in_SD', 'raw', 'n_wrong_mapq60'):,} / {g(m, run, 'in_SD', 'raw', 'n_mapq60'):,}",
                     'MAPQ 60 error, in SD (Phred)': f"{pct(g(m, run, 'in_SD', 'raw', 'w_err_mapq60'))} ({phred(g(m, run, 'in_SD', 'raw', 'w_err_mapq60'))})",
                     'MAPQ≥30 error, in SD': pct(g(m, run, 'in_SD', 'raw', 'mapq30_realised_err'))})
    print('### Table 1. Reported MAPQ versus empirical error\n')
    print(md(pd.DataFrame(rows)))
    # Table 2: recalibrator, in-SD
    rows = []
    for run in RUNS:
        r = {'Run': lab(run)}
        for p, nm in (('raw', 'reported'), ('iso', 'MAPQ-only'), ('gbm', 'GBM'), ('gbmcal', 'GBM+iso')):
            r[f'ECE {nm}'] = '%.3f' % g(m, run, 'in_SD', p, 'ECE')
        e0, e1 = g(m, run, 'in_SD', 'raw', 'ECE'), g(m, run, 'in_SD', 'gbmcal', 'ECE')
        r['ECE change (GBM+iso)'] = '%+.0f%%' % (100 * (e1 / e0 - 1))
        r['Brier reported → GBM+iso'] = '%.3f → %.3f' % (g(m, run, 'in_SD', 'raw', 'Brier'), g(m, run, 'in_SD', 'gbmcal', 'Brier'))
        r['AUPRC(wrong) reported → GBM'] = '%.2f → %.2f' % (g(m, run, 'in_SD', 'raw', 'AUPRC_wrong'), g(m, run, 'in_SD', 'gbm', 'AUPRC_wrong'))
        r['Correct kept at ≤0.1% error: reported / GBM / GBM+iso'] = '%s / %s / %s' % (
            pct(g(m, run, 'in_SD', 'raw', 'frac_correct_retained_err1e-03'), 0), pct(g(m, run, 'in_SD', 'gbm', 'frac_correct_retained_err1e-03'), 0),
            pct(g(m, run, 'in_SD', 'gbmcal', 'frac_correct_retained_err1e-03'), 0))
        rows.append(r)
    print('\n### Table 2. Held-out recalibration inside segmental duplications\n')
    print(md(pd.DataFrame(rows)))
    # Table 3: deployment-style selection (label-free), in SD and genome-wide
    rows = []
    for run in RUNS:
        r = {'Run': lab(run)}
        for st, nm in (('in_SD', 'in SD'), ('all', 'genome-wide')):
            for p, pn in (('raw', 'reported'), ('gbmcal', 'GBM+iso')):
                r[f'{nm}: MAPQ≥30 error / kept, {pn}'] = '%s / %s' % (pct(g(m, run, st, p, 'mapq30_realised_err'), 2), pct(g(m, run, st, p, 'mapq30_frac_correct'), 0))
            r[f'{nm}: realised error at 0.1% target, reported / GBM / GBM+iso'] = '%s / %s / %s' % (
                pct(g(m, run, st, 'raw', 'selfsel_realised_err1e-03'), 2), pct(g(m, run, st, 'gbm', 'selfsel_realised_err1e-03'), 2),
                pct(g(m, run, st, 'gbmcal', 'selfsel_realised_err1e-03'), 2))
        rows.append(r)
    print('\n### Table 3. Selecting reads with the score alone (thresholds never see truth labels)\n')
    print(md(pd.DataFrame(rows)))
    # Table 4: diagnosis of MAPQ-60 errors
    rows = []
    for run in RUNS:
        x = dg[dg.run == run].iloc[0]
        rows.append({'Run': lab(run), 'Wrong at MAPQ 60': int(x.n_wrong_mapq60), 'In SD': pct(x.frac_in_sd, 0),
                     'Read scores better at reported locus': f"{int(x.mapped_better)} ({pct(x.frac_mapped_better, 0)})",
                     'Tie': int(x.tie), 'Truth locus better (search error)': int(x.truth_better),
                     'Read does not align at truth locus': int(x.truth_unaligned), 'Median truth-SD identity': '%.3f' % x.median_truth_sd_identity})
    print('\n### Table 4. Why MAPQ-60 placements are wrong (each read realigned to its truth and reported loci ±2 kb)\n')
    print(md(pd.DataFrame(rows)))
    # Table 5: strict truth sensitivity (long reads)
    rows = []
    for run in [r for r in RUNS if r.split('.')[1] in ('mm2_hifi', 'mm2_lrhq', 'mm2_ont')]:
        rows.append({'Run': lab(run),
                     'MAPQ 60 in-SD error: ≥50% / ≥90% colinear truth': '%s / %s' % (pct(g(m, run, 'in_SD', 'raw', 'w_err_mapq60')), pct(g(ms, run, 'in_SD', 'raw', 'w_err_mapq60'))),
                     'ECE reported / MAPQ-only / GBM / GBM+iso (strict)': ' / '.join('%.3f' % g(ms, run, 'in_SD', p, 'ECE') for p in ('raw', 'iso', 'gbm', 'gbmcal')),
                     'Correct kept at ≤0.1%: reported / GBM / GBM+iso (strict)': ' / '.join(pct(g(ms, run, 'in_SD', p, 'frac_correct_retained_err1e-03'), 0) for p in ('raw', 'gbm', 'gbmcal')),
                     'Realised error at 0.1% target: MAPQ-only / GBM+iso (strict)': '%s / %s' % (pct(g(ms, run, 'in_SD', 'iso', 'selfsel_realised_err1e-03'), 2), pct(g(ms, run, 'in_SD', 'gbmcal', 'selfsel_realised_err1e-03'), 2))})
    print('\n### Table 5. Long reads under a stricter truth (≥90% of the read colinear with its orthologous locus)\n')
    print(md(pd.DataFrame(rows)))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
