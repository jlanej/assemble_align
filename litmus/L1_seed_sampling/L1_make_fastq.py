"""L1_make_fastq.py -- write the simulated read sets as FASTQ with paftools.js-mapeval truth names.

usage: python L1_make_fastq.py DATA_DIR OUT_DIR
Read names: <set>_<i>!chr20!<start>!<end>!<strand>  (0-based start on CHM13 chr20).
CHM13-derived reads: truth = source interval.  NA19240-derived reads: truth = min/max of the lifted
CHM13 coordinates of the read's source bases (reads with <50% liftable bases or with a lifted span
> 1.5x the read length are written to a separate *_unlifted.fq and not scored).
Uses the same seeds as L1_run.py, so the reads are identical to those used for seed metrics.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import L1_schemes as S
import L1_simulate as SIM
import L1_eval as E
import L1_run as R

ACGT = np.frombuffer(b'ACGTN', np.uint8)


def write(rs, name, lift, fo, fu):
    nsc = nun = 0
    for i in range(len(rs['offs'])):
        a, L = rs['offs'][i], rs['lens'][i]
        seq = ACGT[rs['codes'][a:a + L]].tobytes().decode()
        src = rs['src'][a:a + L]
        src = src[src >= 0]
        o = int(rs['meta']['orient'][i])
        if lift is None:
            st, en, strand, ok = int(src.min()), int(src.max()) + 1, '+-'[o], True
        else:
            lp, ls = lift
            t = lp[src]
            good = t >= 0
            ok = good.mean() >= 0.5
            if ok:
                tt = t[good]
                st, en = int(tt.min()), int(tt.max()) + 1
                so = int(np.bincount(ls[src[good]].astype(np.int64), minlength=2).argmax())
                strand = '+-'[o ^ so]
                ok = (en - st) <= 1.5 * L
        if ok:
            fo.write(f'@{name}_{i}!chr20!{st}!{en}!{strand}\n{seq}\n+\n{"I" * L}\n'); nsc += 1
        else:
            fu.write(f'@{name}_{i}\n{seq}\n+\n{"I" * L}\n'); nun += 1
    return nsc, nun


def main(data, out):
    os.makedirs(out, exist_ok=True)
    seqs, lifts = R.load_inputs(data)
    names = ['ill_self', 'ill_pat', 'ill_mat', 'hifi_self', 'hifi_pat', 'hifi_mat', 'ont_self', 'ont_pat', 'ont_mat']
    rsets, cal = R.make_readsets(seqs, names)
    for tech in ('ill', 'hifi', 'ont'):
        with open(os.path.join(out, f'{tech}.fq'), 'w') as fo, open(os.path.join(out, f'{tech}_unlifted.fq'), 'w') as fu:
            for src in ('self', 'pat', 'mat'):
                nm = f'{tech}_{src}'
                n1, n2 = write(rsets[nm], nm, None if src == 'self' else lifts[src], fo, fu)
                print(nm, 'scored', n1, 'unlifted', n2, flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
