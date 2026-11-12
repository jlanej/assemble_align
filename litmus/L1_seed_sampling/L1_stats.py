"""L1_stats.py -- summarise L1 results with paired bootstrap 95% CIs vs the minimap2 baseline.

usage: python L1_stats.py RESULTS_DIR OUT_TSV
Inputs : RESULTS_DIR/index_chm13.tsv, index_grch38.tsv (optional), blocks_chm13.npz (optional),
         perread/<preset>__<scheme>__<readset>.npz
Output : one row per (preset, scheme, read group) with absolute metrics and % change vs 'mm2'
         (bootstrap over reads, B=1000, paired: the same resampled reads for scheme and baseline;
         index metrics: bootstrap over 1-Mb blocks of the reference).
Strata: ':uniq' = reads with <50% baseline seeds above mid_occ; ':rep' = the rest (satellites etc.).
Reads whose truth is undefined (<50% of baseline seeds liftable) are excluded for all schemes.
Read groups: self = CHM13-derived reads (sequencing error only); hap = NA19240 pat+mat reads
(sequencing error + real divergence); asm = error-free 20-kb NA19240 haplotype tiles.
"""
import os, sys, glob
import numpy as np

B = 1000
RNG = np.random.default_rng(2026)


def load_group(res, preset, scheme, sets):
    acc = {}
    for rn in sets:
        f = os.path.join(res, 'perread', f'{preset.replace(":", "_")}__{scheme}__{rn}.npz')
        if not os.path.exists(f):
            return None
        z = np.load(f)
        for k in z.files:
            acc.setdefault(k, []).append(z[k])
    return {k: np.concatenate(v).astype(np.float64) for k, v in acc.items()}


def per_read_stats(d):
    L = d['lens']
    return dict(n_reads=len(L), seeds_per_read=d['n_seeds'].mean(), true_per_read=d['n_true'].mean(),
                true_per_kb=d['n_true'].sum() / L.sum() * 1e3, frac_ge1=(d['n_true'] >= 1).mean(),
                frac_ge3=(d['n_true'] >= 3).mean(), anchors_per_read=d['anchors'].mean(),
                true_frac_of_anchors=d['n_true'].sum() / max(1.0, d['anchors'].sum()),
                hf_frac_seeds=d['n_hf'].sum() / max(1.0, d['n_seeds'].sum()),
                lost_to_hf_per_read=(d['n_true_any'] - d['n_true']).mean(),
                gap_mean=d['max_gap'].mean(), gap_median=np.median(d['max_gap']),
                gap_p95=np.percentile(d['max_gap'], 95))


def boot_ratio(x, y, n, stat='mean'):
    """Paired bootstrap of stat(x)/stat(y)-1 (x = scheme, y = baseline), batched."""
    out = []
    for b0 in range(0, B, 50):
        idx = RNG.integers(0, n, (min(50, B - b0), n))
        if stat == 'mean':
            xs = x[idx].mean(1); ys = y[idx].mean(1)
        else:   # 95th percentile
            xs = np.percentile(x[idx], 95, axis=1); ys = np.percentile(y[idx], 95, axis=1)
        out.append(xs / np.maximum(ys, 1e-12) - 1)
    out = np.concatenate(out)
    return np.percentile(out, [2.5, 97.5])


def boot_diff_frac(x, y, n, thr):
    out = []
    for b0 in range(0, B, 50):
        idx = RNG.integers(0, n, (min(50, B - b0), n))
        out.append((x[idx] >= thr).mean(1) - (y[idx] >= thr).mean(1))
    return np.percentile(np.concatenate(out), [2.5, 97.5])


def read_tsv(path):
    if not os.path.exists(path):
        return []
    lines = open(path).read().strip().split('\n')
    keys = lines[0].split('\t')
    return [dict(zip(keys, l.split('\t'))) for l in lines[1:]]


GROUPS = {'sr': {'self': ['ill_self'], 'hap': ['ill_pat', 'ill_mat']},
          'map-ont': {'self': ['ont_self'], 'hap': ['ont_pat', 'ont_mat'], 'pbself': ['ontpb_self'], 'pbhap': ['ontpb_pat', 'ontpb_mat']},
          'map-hifi': {'self': ['hifi_self'], 'hap': ['hifi_pat', 'hifi_mat'], 'asm': ['asm_pat', 'asm_mat']},
          'asm20': {'asm': ['asm_pat', 'asm_mat']}, 'lr:hqae': {'asm': ['asm_pat', 'asm_mat']}}


def _group_row(row, grp, n_excl, d, d0, sch):
    st = per_read_stats(d); st0 = per_read_stats(d0)
    rr = dict(row); rr['group'] = grp; rr['n_excluded_no_truth'] = n_excl
    for k_, v in st.items():
        rr[k_] = v
    n = len(d['lens'])
    for met, arr in [('true_per_read', 'n_true'), ('anchors_per_read', 'anchors'), ('seeds_per_read', 'n_seeds'),
                     ('gap_mean', 'max_gap')]:
        rr[met + '_pct'] = 100 * (st[met] / max(st0[met], 1e-12) - 1)
        if sch != 'mm2':
            lo, hi = boot_ratio(d[arr], d0[arr], n)
            rr[met + '_pct_ci'] = f'[{100*lo:.2f}, {100*hi:.2f}]'
    rr['gap_p95_pct'] = 100 * (st['gap_p95'] / max(st0['gap_p95'], 1e-12) - 1)
    for thr, nm in [(1, 'frac_ge1'), (3, 'frac_ge3')]:
        rr[nm + '_diff_pp'] = 100 * (st[nm] - st0[nm])
        if sch != 'mm2':
            lo, hi = boot_diff_frac(d['n_true'], d0['n_true'], n, thr)
            rr[nm + '_diff_pp_ci'] = f'[{100*lo:.3f}, {100*hi:.3f}]'
    return rr


def main(res, out):
    idx = read_tsv(os.path.join(res, 'index_chm13.tsv'))
    idx38 = {(r['preset'], r['scheme']): r for r in read_tsv(os.path.join(res, 'index_grch38.tsv'))}
    blocks = np.load(os.path.join(res, 'blocks_chm13.npz')) if os.path.exists(os.path.join(res, 'blocks_chm13.npz')) else None
    rows = []
    for r in idx:
        preset, sch = r['preset'], r['scheme']
        base = next(x for x in idx if x['preset'] == preset and x['scheme'] == 'mm2')
        row = dict(preset=preset, scheme=sch, k=r['k'], w=r['w'], s=r.get('s', ''),
                   density=float(r['density']), dfactor=float(r['dfactor']), entries=int(r['entries']),
                   entries_pct=100 * (int(r['entries']) / int(base['entries']) - 1),
                   mid_occ=int(r['mid_occ']), frac_entries_hf=float(r['frac_entries_hf']),
                   ref_max_gap=int(r['ref_max_gap']))
        key = f'{preset}__{sch}'
        if blocks is not None and key in blocks.files and f'{preset}__mm2' in blocks.files:
            xb, yb = blocks[key].astype(float), blocks[f'{preset}__mm2'].astype(float)
            m = yb > 0
            lo, hi = boot_ratio(xb[m], yb[m], int(m.sum()))
            row['entries_pct_ci'] = f'[{100*lo:.2f}, {100*hi:.2f}]'
        if (preset, sch) in idx38:
            g = idx38[(preset, sch)]; gb = idx38[(preset, 'mm2')]
            row['grch38_density'] = float(g['density'])
            row['grch38_entries_pct'] = 100 * (int(g['entries']) / int(gb['entries']) - 1)
            row['grch38_frac_entries_hf'] = float(g['frac_entries_hf'])
        for grp, sets in GROUPS[preset].items():
            d = load_group(res, preset, sch, sets); d0 = load_group(res, preset, 'mm2', sets)
            if d is None or d0 is None:
                continue
            # truth-defined reads: >=50% of the baseline's seeds have a defined (lifted) expected position
            mask = (d0['n_seeds'] > 0) & (d0['liftable'] >= 0.5 * d0['n_seeds'])
            n_excl = int((~mask).sum())
            dA = {k_: v[mask] for k_, v in d.items()}; d0A = {k_: v[mask] for k_, v in d0.items()}
            rep = d0A['n_hf'] >= 0.5 * d0A['n_seeds']     # repeat stratum, defined on the baseline
            for sub, m in [('', None), (':uniq', ~rep), (':rep', rep)]:
                if m is None:
                    d, d0 = dA, d0A
                else:
                    if m.sum() < 5:
                        continue
                    d = {k_: v[m] for k_, v in dA.items()}; d0 = {k_: v[m] for k_, v in d0A.items()}
                rows.append(_group_row(row, grp + sub, n_excl, d, d0, sch))
    keys = []
    for r in rows:
        for k_ in r:
            if k_ not in keys:
                keys.append(k_)
    with open(out, 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(k_) is None else (f'{r[k_]:.6g}' if isinstance(r.get(k_), (float, np.floating)) else str(r[k_])) for k_ in keys) + '\n')
    print('wrote', out, len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
