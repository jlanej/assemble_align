"""L1_validate.py -- correctness checks for L1_schemes.py.

1. chunked == unchunked sampling (all schemes) on a real-sequence slice.
2. strand symmetry: sample(seq) vs sample(revcomp(seq)) mapped back (fraction of positions shared).
3. density on an i.i.d. uniform random sequence vs. theory / published values:
     random minimizer 2/(w+1); mod-minimizer (2+(k-t)/w)/(w+k-t+1) [Groot Koerkamp & Pibiri 2024];
     miniception / open-closed at (w=5,k=11,s=6) = 0.2929 / 0.2864 [Groot Koerkamp, Liu & Pibiri 2025];
     lower bound ceil((w+k)/w)/(w+k) [Kille et al. 2024].
Usage: python L1_validate.py <chr20.fa> <out.tsv> [harness_stock harness_modoc]
(The minimap2 bit-identity check uses harness/sketch_harness.c; see L1_run_harness.sh.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys, time
import numpy as np
import L1_schemes as S


def sym_fraction(codes, scheme, k, w, s=None, greedy=None):
    n = len(codes)
    p1, _, c1 = S.sample(codes, scheme, k, w, s=s, greedy=greedy)
    rc = S.revcomp_codes(codes)
    p2, _, c2 = S.sample(rc, scheme, k, w, s=s, greedy=greedy)
    p2m = np.sort(n + k - 2 - p2)
    inter = np.intersect1d(p1, p2m).size
    return inter / max(1, len(p1)), inter / max(1, len(p2m)), len(p1), len(p2m)


def main(fa, out):
    rows = []
    name, seq = S.read_fasta(fa)[0]
    codes = S.encode(seq)
    go15, kp15 = S.load_greedy_order(os.path.join(os.environ.get('L1_GM_DIR', 'src/gm'), 'k15-w10.txt'))
    # --- 1. chunked == unchunked, 3 Mb slice, small chunk to force many boundaries
    sl = codes[30_000_000:33_000_000]
    for sch, k, w, s, g in [('mm2', 15, 10, None, None), ('oc', 19, 19, 4, None), ('ddecyc', 15, 10, None, None),
                            ('mod', 21, 11, None, None), ('modoc', 15, 10, None, None), ('osync', 15, 10, 10, None),
                            ('miniception', 15, 10, None, None), ('greedy', 15, 10, None, (go15, kp15))]:
        a = S.sample(sl, sch, k, w, s=s, greedy=g)
        b = S.sample_chunked(sl, sch, k, w, s=s, greedy=g, chunk=400_000)
        ok = all(np.array_equal(x, y) for x, y in zip(a, b))
        rows.append(('chunk_equal', sch, k, w, s, int(ok), len(a[0])))
        print('chunk_equal', sch, k, w, ok, len(a[0]), len(b[0]))
    # --- 2. strand symmetry on a 2 Mb real slice
    sl = codes[10_000_000:12_000_000]
    for sch, k, w, s, g in [('mm2', 15, 10, None, None), ('random', 15, 10, None, None), ('miniception', 15, 10, None, None),
                            ('oc', 19, 19, 4, None), ('ddecyc', 15, 10, None, None), ('greedy', 15, 10, None, (go15, kp15)),
                            ('mod', 15, 10, None, None), ('mod', 21, 11, None, None), ('modoc', 15, 10, None, None),
                            ('osync', 15, 10, 10, None), ('csync', 15, 10, 5, None)]:
        f1, f2, n1, n2 = sym_fraction(sl, sch, k, w, s, g)
        rows.append(('strand_sym', sch, k, w, s, round(min(f1, f2), 6), n1))
        print('strand_sym', sch, k, w, s, round(f1, 6), round(f2, 6), n1, n2)
    # --- 3. densities on random sequence (10 Mb, seed 1), forward (paper semantics) and canonical
    rng = np.random.default_rng(1)
    rnd = rng.integers(0, 4, 10_000_000).astype(np.uint8)
    nk = lambda k: len(rnd) - k + 1
    tests = [('random', 5, 11, None), ('miniception', 5, 11, 6), ('oc', 5, 11, 6),
             ('random', 10, 15, None), ('mod', 10, 15, None), ('random', 11, 21, None), ('mod', 11, 21, None),
             ('random', 19, 19, None), ('oc', 19, 19, 4), ('miniception', 19, 19, None), ('ddecyc', 19, 19, None),
             ('ddecyc', 10, 15, None), ('modoc', 10, 15, None), ('greedy', 10, 15, None)]
    for sch, w, k, s in tests:
        g = (go15, kp15) if sch == 'greedy' else None
        if sch == 'miniception' and s is not None:
            # miniception with explicit r=s (published check value): closed syncmers first, forward semantics
            so, sc = S.smer_classes(rnd, k, s, False, 0x5bd1e995)
            fwd, rev, valid = S.kmer_codes(rnd, k)
            fh = S.fwd_hash(fwd, valid, k)
            rank = S.class_rank(fh, valid, np.where(sc == 1, 0, 1).astype(np.uint8))
            pf = S.window_min_select_fwd(rank, valid, w, k, len(rnd))
            pc = None
        else:
            pf = S.sample_forward(rnd, sch, k, w, s=s, greedy=g)
            pc = S.sample(rnd, sch, k, w, s=s, greedy=g)[0]
        dfwd = len(pf) / nk(k)
        dcan = len(pc) / nk(k) if pc is not None else float('nan')
        t = S.mod_t(k, w)
        theory = {'random': 2 / (w + 1), 'mod': (2 + (k - t) // w) / (w + k - t + 1)}.get(sch, float('nan'))
        if (sch, w, k) == ('miniception', 5, 11): theory = 0.2929
        if (sch, w, k) == ('oc', 5, 11): theory = 0.2864
        lb = int(np.ceil((w + k) / w)) / (w + k)
        rows.append(('density_random', sch, k, w, s, round(dfwd, 5), round(dcan, 5), round(theory, 5), round(lb, 5)))
        print('density', sch, 'w', w, 'k', k, 's', s, 'fwd', round(dfwd, 5), 'canon', round(dcan, 5),
              'theory/published', round(theory, 5), 'LB', round(lb, 5))
    with open(out, 'w') as f:
        f.write('test\tscheme\tk\tw\ts\tvalue1\tvalue2\tvalue3\tvalue4\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


def validate_harness(fa, harness_stock, harness_modoc, out):
    """Bit-identity of the numba port vs minimap2's own sketch.c (stock) and of the Python mod-OC vs the
    patched sketch.c (harness built from build/mm2-modoc), on every sequence of `fa` (first record)."""
    import subprocess, tempfile, os
    name, seq = S.read_fasta(fa)[0]
    codes = S.encode(seq)
    rows = []
    with tempfile.TemporaryDirectory(dir='.') as d:
        for k, w in [(21, 11), (15, 10), (19, 19), (19, 10), (25, 51), (15, 5)]:
            b = os.path.join(d, 'x.bin')
            subprocess.run([harness_stock, fa, str(k), str(w), b], check=True, stderr=subprocess.DEVNULL)
            ref = np.fromfile(b, dtype=np.uint64).reshape(-1, 2)
            x, y = S.mm_sketch_generic(codes, w, k, np.empty(1, np.uint64), False)
            same = len(x) == len(ref) and np.array_equal(x, ref[:, 0]) and np.array_equal(y, ref[:, 1])
            rows.append(('mm2_port_vs_sketch.c', 'mm2', k, w, '', int(same), len(ref)))
            if harness_modoc and (k, w) != (15, 5):
                subprocess.run([harness_modoc, fa, str(k), str(w), b, '8'], check=True, stderr=subprocess.DEVNULL)
                c = np.fromfile(b, dtype=np.uint64).reshape(-1, 2)
                cpos = np.sort((c[:, 1] >> np.uint64(1)).astype(np.int64))
                ppos = S.sample2_chunked(codes, 'modoc', k, w, chunk=3_000_000)[0]
                rows.append(('modoc_python_vs_patched_sketch.c', 'modoc', k, w, '', int(np.array_equal(cpos, np.sort(ppos))), len(cpos)))
            print(rows[-1], flush=True)
    # hash quality on i.i.d. random sequence (4 Mb): hash64 vs width-adapted invertible mixer vs splitmix56
    rng = np.random.default_rng(1)
    rnd = rng.integers(0, 4, 4_000_000).astype(np.uint8)
    for (w, k) in [(5, 11), (10, 15), (5, 15), (11, 21), (19, 19)]:
        fwd, rev, valid = S.kmer_codes(rnd, k)
        n = len(rnd) - k + 1
        x1, _ = S.mm_sketch_generic(rnd, w, k, np.empty(1, np.uint64), False)
        x2, _ = S.mm_sketch_generic(rnd, w, k, S.canon_hash_w(fwd, rev, valid, k), True)
        x3, _ = S.mm_sketch_generic(rnd, w, k, S.canon_mix(fwd, rev, valid), True)
        rows.append(('hash_density_random', 'hash64|hash_w|splitmix56', k, w, '', f'{len(x1)/n:.5f}|{len(x2)/n:.5f}|{len(x3)/n:.5f}', f'{2/(w+1):.5f}'))
        print(rows[-1], flush=True)
    # open-closed minimizer: forward vs canonical density as a function of s (k-s parity), 3 Mb random
    sub = rnd[:3_000_000]
    for (w, k) in [(19, 19), (51, 25), (11, 21), (10, 15), (10, 19)]:
        fwd, rev, valid = S.kmer_codes(sub, k)
        chf = S.fwd_mix(fwd, valid); chc = S.canon_mix(fwd, rev, valid)
        for s in range(3, min(k, 14)):
            vals = []
            for canonical in (False, True):
                so, sc = S.smer_classes_mix(sub, k, s, canonical)
                cls = np.where(so == 1, 0, np.where(sc == 1, 1, 2)).astype(np.uint8)
                rank = S.class_rank(chc if canonical else chf, valid, cls)
                if canonical:
                    x, _ = S.mm_sketch_generic(sub, w, k, rank, True); n_s = len(x)
                else:
                    n_s = len(S.window_min_select_fwd(rank, valid, w, k, len(sub)))
                vals.append(n_s / (len(sub) - k + 1))
            rows.append(('oc_density_vs_s', 'oc fwd|canonical', k, w, s, f'{vals[0]:.5f}|{vals[1]:.5f}', 'k-s even' if (k - s) % 2 == 0 else 'k-s odd'))
        print('oc scan', w, k, flush=True)
    rows.append(('hash_w_bijection_k8_to_12', 'hash_w', '8-12', '', '', int(all(S.hash_w_is_bijection(k) for k in (8, 10, 11, 12))), ''))
    with open(out, 'a') as f:
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
    if len(sys.argv) > 3:   # optional: paths to harness binaries (stock, modoc)
        validate_harness(sys.argv[1], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None, sys.argv[2])
