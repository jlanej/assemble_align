"""L1_tables.py -- render the markdown tables of L1_sampling_schemes.md from the result TSVs.
usage: python L1_tables.py L1_results.tsv RESULTS_DIR OUT.md
"""
import os, sys
import pandas as pd

ORDER = ['mm2', 'mm2_1pick', 'hashw', 'rand_mix', 'miniception', 'ddecyc', 'greedy', 'mod', 'modoc', 'oc', 'osync', 'csync']
LAB = {'mm2': 'minimap2 2.31 (baseline)', 'mm2_1pick': 'same order, 1 k-mer/window', 'hashw': 'minimap2 order, width-adapted invertible hash', 'rand_mix': 'random, strong hash',
       'miniception': 'miniception', 'ddecyc': 'double decycling', 'greedy': 'GreedyMini (min over strands)',
       'mod': 'mod-minimizer', 'modoc': 'open-closed mod-minimizer', 'oc': 'open-closed minimizer',
       'osync': 'open syncmer', 'csync': 'closed syncmer'}
PRES = ['sr', 'map-ont', 'map-hifi', 'asm20', 'lr:hqae']


def f(x, nd=2):
    return '' if pd.isna(x) else f'{x:.{nd}f}'


def main(res_tsv, resdir, out):
    R = pd.read_csv(res_tsv, sep='\t')
    L = []
    # Table 1: index
    L.append('| preset (k,w) | scheme | s | density | d·(w+1) | entries | Δentries vs mm2 (95% CI) | mid_occ | high-freq entries | GRCh38 Δentries |')
    L.append('|---|---|---|---|---|---|---|---|---|---|')
    ix = R.drop_duplicates(['preset', 'scheme'])
    for p in PRES:
        sub = ix[ix.preset == p].copy()
        sub['o'] = sub.scheme.map({s: i for i, s in enumerate(ORDER)})
        for _, r in sub.sort_values('o').iterrows():
            ci = '' if r.scheme == 'mm2' else f" {r.entries_pct_ci}"
            g38 = '' if pd.isna(r.get('grch38_entries_pct', float('nan'))) else f"{r.grch38_entries_pct:+.1f}%"
            L.append(f"| {p} ({int(r.k)},{int(r.w)}) | {LAB[r.scheme]} | {'' if pd.isna(r.s) else int(r.s)} | {r.density:.4f} | {r.dfactor:.3f} | "
                     f"{int(r.entries):,} | {r.entries_pct:+.1f}%{ci} | {int(r.mid_occ)} | {100*r.frac_entries_hf:.1f}% | {g38} |")
    L.append('')
    # Table 2: per-read, main groups
    for grp_title, grps in [('NA19240-derived reads and haplotype tiles (sequencing error + real divergence)', ['hap', 'asm', 'pbhap']),
                            ('CHM13-derived reads (sequencing error only)', ['self', 'pbself'])]:
        L.append(f'**{grp_title}.**')
        L.append('')
        L.append('| preset | reads | scheme | true seeds/read | Δ true (95% CI) | Δ anchors (95% CI) | ≥1 true | ≥3 true | Δ≥3 pp (95% CI) | longest gap mean / p95 (bp) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|')
        for p in PRES:
            for g in grps:
                sub = R[(R.preset == p) & (R.group == g)].copy()
                if sub.empty:
                    continue
                sub['o'] = sub.scheme.map({s: i for i, s in enumerate(ORDER)})
                for _, r in sub.sort_values('o').iterrows():
                    if r.scheme == 'mm2_1pick':
                        continue
                    dt = '—' if r.scheme == 'mm2' else f"{r.true_per_read_pct:+.1f}% {r.true_per_read_pct_ci}"
                    da = '—' if r.scheme == 'mm2' else f"{r.anchors_per_read_pct:+.1f}% {r.anchors_per_read_pct_ci}"
                    d3 = '—' if r.scheme == 'mm2' else f"{r.frac_ge3_diff_pp:+.2f} {r.frac_ge3_diff_pp_ci}"
                    L.append(f"| {p} | {g} (n={int(r.n_reads)}) | {LAB[r.scheme]} | {r.true_per_read:.1f} | {dt} | {da} | {100*r.frac_ge1:.2f}% | "
                             f"{100*r.frac_ge3:.2f}% | {d3} | {r.gap_mean:.0f} / {r.gap_p95:.0f} |")
        L.append('')
    # Table 3: strata for the focal schemes
    L.append('**Unique vs repeat strata (NA19240 reads/tiles; stratum defined on the baseline: ≥50% of seeds above mid_occ = repeat).**')
    L.append('')
    L.append('| preset | stratum | scheme | n | Δ true (95% CI) | ≥1 true | Δ≥1 pp (95% CI) | Δ≥3 pp (95% CI) | longest gap p95 (bp) |')
    L.append('|---|---|---|---|---|---|---|---|---|')
    for p in PRES:
        for g in ['hap:uniq', 'hap:rep', 'asm:uniq', 'asm:rep']:
            for sch in ['mm2', 'rand_mix', 'mod', 'modoc', 'oc', 'csync', 'osync']:
                sub = R[(R.preset == p) & (R.group == g) & (R.scheme == sch)]
                if sub.empty:
                    continue
                r = sub.iloc[0]
                dt = '—' if sch == 'mm2' else f"{r.true_per_read_pct:+.1f}% {r.true_per_read_pct_ci}"
                d1 = '—' if sch == 'mm2' else f"{r.frac_ge1_diff_pp:+.2f} {r.frac_ge1_diff_pp_ci}"
                d3 = '—' if sch == 'mm2' else f"{r.frac_ge3_diff_pp:+.2f} {r.frac_ge3_diff_pp_ci}"
                L.append(f"| {p} | {g} | {LAB[sch]} | {int(r.n_reads)} | {dt} | {100*r.frac_ge1:.2f}% | {d1} | {d3} | {r.gap_p95:.0f} |")
    L.append('')
    # Table 4: patch evaluation
    pe = os.path.join(resdir, 'patch', 'patch_eval.tsv')
    if os.path.exists(pe):
        P = pd.read_csv(pe, sep='\t')
        L.append('| preset | variant | reads | .mmi bytes | index CPU s | index RSS GB | map CPU s | map RSS GB | MAPQ60 mapped (err) | MAPQ≥1 mapped (err) | all mapped (err) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|')
        for _, r in P.iterrows():
            L.append(f"| {r.preset} | {'stock 2.31' if r.variant == 'stock' else '2.31 + --mod-oc'} | {int(r.reads)} | {int(r.mmi_bytes):,} | {r.idx_cpu_s:.2f} | {r.idx_rss_gb:.3f} | "
                     f"{r.map_cpu_s:.2f} | {r.map_rss_gb:.3f} | {int(r.mapped_top)} ({100*r.err_top:.3f}%) | {int(r.mapped_q1)} ({100*r.err_q1:.3f}%) | {int(r.mapped_all)} ({100*r.err_all:.3f}%) |")
        L.append('')
    open(out, 'w').write('\n'.join(L) + '\n')
    print('wrote', out, len(L), 'lines')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
