"""L1_run.py -- litmus L1 driver: sampling schemes x minimap2 presets on chr20.

usage: python L1_run.py DATA_DIR OUT_DIR PART
  PART = main    : CHM13 chr20 indexes for every preset x scheme + per-read seed metrics for the
                   simulated read sets (self = CHM13-derived, hap = NA19240 pat/mat-derived) and
                   NA19240 haplotype tiles (asm presets).  Writes OUT/index_chm13.tsv, OUT/perread/*.npz
         grch38  : GRCh38 chr20 index statistics for every preset x scheme (OUT/index_grch38.tsv)
         curves  : equal-density curves (vary w / s) on hap read sets (OUT/curves.tsv)
DATA_DIR must contain chm13_chr20.fa, hg38_chr20.fa, NA19240_{pat,mat}_chr20.fa and the PAFs
NA19240_{pat,mat}_vs_chm13_chr20.paf (minimap2 -x asm5 -c), and src/gm/*.txt GreedyMini orders.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os, sys, time, json
import numpy as np
import L1_schemes as S
import L1_simulate as SIM
import L1_eval as E

GM_DIR = os.environ.get('L1_GM_DIR', 'src/gm')
DATA_DIR = 'data'
GREEDY_FOR = {'sr': 'k16-w11', 'map-ont': 'k15-w10', 'map-hifi': 'k15-w19', 'asm20': 'k15-w10'}
OSYNC_S = {'sr': 15, 'map-ont': 11, 'map-hifi': 11, 'asm20': 15}
CSYNC_S = {'sr': 10, 'map-ont': 5, 'asm20': 9}
PRESETS = ['sr', 'map-ont', 'map-hifi', 'asm20', 'lr:hqae']
READSETS = {'sr': ['ill_self', 'ill_pat', 'ill_mat'], 'map-ont': ['ont_self', 'ont_pat', 'ont_mat', 'ontpb_self', 'ontpb_pat', 'ontpb_mat'],
            'map-hifi': ['hifi_self', 'hifi_pat', 'hifi_mat', 'asm_pat', 'asm_mat'],
            'asm20': ['asm_pat', 'asm_mat'], 'lr:hqae': ['asm_pat', 'asm_mat']}
SEEDS = {'ill_self': 101, 'ill_pat': 102, 'ill_mat': 103, 'hifi_self': 201, 'hifi_pat': 202, 'hifi_mat': 203,
         'ont_self': 301, 'ont_pat': 302, 'ont_mat': 303}
NREADS = {'ill_self': 20000, 'ill_pat': 10000, 'ill_mat': 10000, 'hifi_self': 1500, 'hifi_pat': 750,
          'hifi_mat': 750, 'ont_self': 800, 'ont_pat': 400, 'ont_mat': 400}


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def scheme_list(preset):
    L = [('mm2', None), ('mm2_1pick', None), ('rand_mix', None), ('miniception', None), ('ddecyc', None)]
    if preset in GREEDY_FOR:
        L.append(('greedy', None))
    L += [('mod', None), ('modoc', None), ('oc', None)]
    if preset in OSYNC_S:
        L.append(('osync', OSYNC_S[preset]))
    if preset in CSYNC_S:
        L.append(('csync', CSYNC_S[preset]))
    return L


def load_inputs(data):
    seqs = {}
    for nm, f in [('chm13', 'chm13_chr20.fa'), ('pat', 'NA19240_pat_chr20.fa'), ('mat', 'NA19240_mat_chr20.fa')]:
        seqs[nm] = S.encode(S.read_fasta(os.path.join(data, f))[0][1])
    lifts = {h: E.load_lift(os.path.join(data, f'NA19240_{h}_vs_chm13_chr20.paf'), len(seqs[h])) for h in ('pat', 'mat')}
    return seqs, lifts


def make_readsets(seqs, which):
    c_hifi = SIM.calibrate_hp(seqs['chm13'], SIM.MODELS['hifi']['hp'])
    c_ont = SIM.calibrate_hp(seqs['chm13'], SIM.MODELS['ont']['hp'])
    rsets = {}
    for nm in which:
        kind, src = nm.split('_')
        sc = seqs['chm13' if src == 'self' else src]
        sid = {'self': 0, 'pat': 1, 'mat': 2}[src]
        if kind == 'asm':
            rsets[nm] = SIM.tiles(sc, 20000, src_id=sid, step=40000)
        elif kind == 'ontpb':
            fa = os.path.join(DATA_DIR, {'self': 'chm13_chr20.fa', 'pat': 'NA19240_pat_chr20.fa', 'mat': 'NA19240_mat_chr20.fa'}[src])
            rsets[nm] = SIM.simulate_pbsim3(fa, sc, 0.36 if src == 'self' else 0.18, {'self': 401, 'pat': 402, 'mat': 403}[src], src_id=sid)
        else:
            k2 = {'ill': 'illumina', 'hifi': 'hifi', 'ont': 'ont'}[kind]
            rsets[nm] = SIM.simulate(sc, k2, NREADS[nm], SEEDS[nm], c_hp={'hifi': c_hifi, 'ont': c_ont}.get(k2), src_id=sid)
    return rsets, dict(c_hp_hifi=c_hifi, c_hp_ont=c_ont)


def greedy_for(preset):
    if preset not in GREEDY_FOR:
        return None
    return S.load_greedy_order(os.path.join(GM_DIR, GREEDY_FOR[preset] + '.txt'))


def run_main(data, out):
    os.makedirs(os.path.join(out, 'perread'), exist_ok=True)
    seqs, lifts = load_inputs(data)
    allrs = sorted({r for p in PRESETS for r in READSETS[p]})
    rsets, cal = make_readsets(seqs, allrs)
    json.dump(dict(cal, readsets={k: dict(n=int(len(v['offs'])), bases=int(v['lens'].sum())) for k, v in rsets.items()}),
              open(os.path.join(out, 'readsets.json'), 'w'), indent=1)
    log('read sets ready', {k: len(v['offs']) for k, v in rsets.items()})
    rows = []
    for preset in PRESETS:
        k, w = E.PRESET_KW[preset]
        g = greedy_for(preset)
        for sch, s in scheme_list(preset):
            t0 = time.time()
            ix = E.Index(seqs['chm13'], sch, k, w, preset, s=s, greedy=g)
            row = ix.summary(); row['s'] = s; row['ref'] = 'chm13_chr20'; row['t_index_s'] = round(time.time() - t0, 1)
            rows.append(row)
            for rn in READSETS[preset]:
                rs = rsets[rn]
                src = rn.split('_')[1]
                lift = None if src == 'self' else lifts[src]
                res = E.eval_reads(ix, rs, sch, k, w, lift=lift, s=s, greedy=g)
                np.savez_compressed(os.path.join(out, 'perread', f'{preset.replace(":", "_")}__{sch}__{rn}.npz'),
                                    lens=rs['lens'], **res)
            log(preset, sch, 'density', round(ix.density, 5), 'entries', ix.n_entries, 'mid_occ', ix.mid_occ,
                'hf', round(ix.frac_entries_hf, 4), f'{time.time()-t0:.1f}s')
            del ix
            write_tsv(rows, os.path.join(out, 'index_chm13.tsv'))


def run_grch38(data, out):
    os.makedirs(out, exist_ok=True)
    hg = S.encode(S.read_fasta(os.path.join(data, 'hg38_chr20.fa'))[0][1])
    rows = []
    for preset in PRESETS:
        k, w = E.PRESET_KW[preset]
        g = greedy_for(preset)
        for sch, s in scheme_list(preset):
            t0 = time.time()
            ix = E.Index(hg, sch, k, w, preset, s=s, greedy=g)
            row = ix.summary(); row['s'] = s; row['ref'] = 'grch38_chr20'; row['t_index_s'] = round(time.time() - t0, 1)
            rows.append(row)
            log('grch38', preset, sch, round(ix.density, 5), ix.mid_occ, round(ix.frac_entries_hf, 4))
            del ix
            write_tsv(rows, os.path.join(out, 'index_grch38.tsv'))


CURVE_W = {'sr': [6, 8, 11, 14, 18], 'map-ont': [6, 8, 10, 13, 16], 'map-hifi': [10, 14, 19, 25, 32]}
CURVE_RS = {'sr': ['ill_pat', 'ill_mat'], 'map-ont': ['ont_pat', 'ont_mat'], 'map-hifi': ['hifi_pat', 'hifi_mat']}


def run_curves(data, out):
    os.makedirs(out, exist_ok=True)
    seqs, lifts = load_inputs(data)
    rsets, cal = make_readsets(seqs, sorted({r for v in CURVE_RS.values() for r in v}))
    rows = []
    for preset, ws in CURVE_W.items():
        k, _ = E.PRESET_KW[preset]
        cfgs = [(sch, w, None) for sch in ('mm2', 'mod', 'modoc', 'oc', 'miniception') for w in ws]
        cfgs += [('osync', 0, s) for s in range(k - 2, 3, -2)][:5]
        cfgs += [('csync', 0, s) for s in range(k - max(ws), k - min(ws) + 1, 2) if s >= 3][:5]
        for sch, w, s in cfgs:
            t0 = time.time()
            wi = w if w else 10
            ix = E.Index(seqs['chm13'], sch, k, wi, preset, s=s)
            agg = dict(preset=preset, scheme=sch, k=k, w=w, s=s, density=ix.density, entries=ix.n_entries,
                       mid_occ=ix.mid_occ, frac_entries_hf=ix.frac_entries_hf)
            acc = {}
            for rn in CURVE_RS[preset]:
                res = E.eval_reads(ix, rsets[rn], sch, k, wi, lift=lifts[rn.split('_')[1]], s=s)
                res['lens'] = rsets[rn]['lens']
                for kk, vv in res.items():
                    acc.setdefault(kk, []).append(vv)
            acc = {kk: np.concatenate(vv) for kk, vv in acc.items()}
            L = acc['lens'].astype(float)
            agg.update(true_per_kb=float(acc['n_true'].sum() / L.sum() * 1000),
                       seeds_per_kb=float(acc['n_seeds'].sum() / L.sum() * 1000),
                       anchors_per_kb=float(acc['anchors'].sum() / L.sum() * 1000),
                       frac_ge1=float((acc['n_true'] >= 1).mean()), frac_ge3=float((acc['n_true'] >= 3).mean()),
                       gap_median=float(np.median(acc['max_gap'])), gap_p95=float(np.percentile(acc['max_gap'], 95)),
                       t_s=round(time.time() - t0, 1))
            rows.append(agg)
            log('curve', preset, sch, w, s, round(ix.density, 5), round(agg['true_per_kb'], 2))
            del ix
            write_tsv(rows, os.path.join(out, 'curves.tsv'))


def run_blocks(data, out):
    """Per-1-Mb counts of sampled positions on CHM13 chr20 (for the index-level block bootstrap)."""
    os.makedirs(out, exist_ok=True)
    chm = S.encode(S.read_fasta(os.path.join(data, 'chm13_chr20.fa'))[0][1])
    nb = (len(chm) + 999_999) // 1_000_000
    res = {}
    for preset in PRESETS:
        k, w = E.PRESET_KW[preset]
        g = greedy_for(preset)
        for sch, s in scheme_list(preset):
            pos, st, cn = S.sample2_chunked(chm, sch, k, w, s=s, greedy=g)
            res[f'{preset}__{sch}'] = np.bincount(pos // 1_000_000, minlength=nb)
            log('blocks', preset, sch, len(pos))
    np.savez_compressed(os.path.join(out, 'blocks_chm13.npz'), **res)


def run_hashw(data, out):
    """L1-B: minimap2 window logic with the width-adapted invertible hash at the k=15 preset (map-ont)."""
    seqs, lifts = load_inputs(data)
    preset = 'map-ont'
    k, w = E.PRESET_KW[preset]
    rsets, cal = make_readsets(seqs, READSETS[preset])
    t0 = time.time()
    ix = E.Index(seqs['chm13'], 'hashw', k, w, preset)
    row = ix.summary(); row['s'] = None; row['ref'] = 'chm13_chr20'; row['t_index_s'] = round(time.time() - t0, 1)
    for rn in READSETS[preset]:
        src = rn.split('_')[1]
        res = E.eval_reads(ix, rsets[rn], 'hashw', k, w, lift=None if src == 'self' else lifts[src])
        np.savez_compressed(os.path.join(out, 'perread', f'{preset}__hashw__{rn}.npz'), lens=rsets[rn]['lens'], **res)
    path = os.path.join(out, 'index_chm13.tsv')
    lines = open(path).read().strip().split('\n'); keys = lines[0].split('\t')
    rows = [dict(zip(keys, l.split('\t'))) for l in lines[1:] if l.split('\t')[keys.index('scheme')] != 'hashw']
    rows.append(row)
    write_tsv(rows, path)
    bp = os.path.join(out, 'blocks_chm13.npz')
    blocks = dict(np.load(bp)) if os.path.exists(bp) else {}
    pos = S.sample2_chunked(seqs['chm13'], 'hashw', k, w)[0]
    blocks[f'{preset}__hashw'] = np.bincount(pos // 1_000_000, minlength=(len(seqs['chm13']) + 999_999) // 1_000_000)
    np.savez_compressed(bp, **blocks)
    log('hashw', preset, 'density', round(ix.density, 5), 'entries', ix.n_entries, 'mid_occ', ix.mid_occ)


def write_tsv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(k) is None else (f'{r[k]:.6g}' if isinstance(r.get(k), float) else str(r[k])) for k in keys) + '\n')


if __name__ == '__main__':
    data, out, part = sys.argv[1], sys.argv[2], sys.argv[3]
    DATA_DIR = data
    {'main': run_main, 'grch38': run_grch38, 'curves': run_curves, 'blocks': run_blocks, 'hashw': run_hashw}[part](data, out)
    import resource
    log('done; peak RSS (GB):', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2))
