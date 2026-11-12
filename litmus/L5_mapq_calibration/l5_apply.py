#!/usr/bin/env python3
"""l5_apply.py - rewrite MAPQ with a trained L5 recalibrator (model_*.pkl from l5_recal.py --save-models).

Input must be grouped by read name (raw aligner output, or `samtools sort -n`). For each read the
evaluated record (mapeval mode 0) and the primary record get
    MAPQ  <- min(60, round(-10 log10 P(wrong)))   (recalibrated; GBM + isotonic post-calibration by default)
    om:i  <- original MAPQ ;  pw:f <- recalibrated P(wrong)
Secondary records are written unchanged. Features are computed exactly as for training.

usage: l5_apply.py --model model_<run>.pkl --sd SD_FILE --ref grch38|chm13 --in ALN.sam|bam --out OUT.bam
"""
import argparse, os, pickle, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # PYTHONSAFEPATH-safe
import numpy as np
import pysam
from l5_common import SDIndex, read_sd_table
from l5_features import new_state, update, rows_from, add_mate_and_sd, mate_of
from l5_recal import add_derived, to_mapq, EPS


def process(batch, model, feats, sd, out):
    states = {}
    for qn, recs in batch:
        for i, r in enumerate(recs):
            k = (qn, mate_of(r))
            if k not in states:
                states[k] = new_state()
            update(states[k], r, rec_id=i)
    df = add_derived(add_mate_and_sd(rows_from(states), sd))
    for c in feats:
        if c not in df.columns:
            df[c] = np.nan
    p = np.full(len(df), np.nan)
    mp = df.mapped.values == 1
    if mp.any():
        if isinstance(model, tuple):            # (GBM fitted on 70%, isotonic calibrator fitted on 30%)
            p[mp] = np.clip(model[1].predict(model[0].predict_proba(df.loc[mp, feats])[:, 1]), EPS, 1 - EPS)
        else:
            p[mp] = np.clip(model.predict_proba(df.loc[mp, feats])[:, 1], EPS, 1 - EPS)
    pw = {(q, m): v for q, m, v in zip(df.qname.values, df.mate.values, p)}
    n_changed = 0
    for qn, recs in batch:
        for i, r in enumerate(recs):
            k = (qn, mate_of(r))
            st = states[k]
            v = pw.get(k, np.nan)
            if (not np.isnan(v)) and (not r.is_unmapped) and (not r.is_secondary) and \
                    (i == st['best_id'] or not r.is_supplementary):
                r.set_tag('om', r.mapping_quality, 'i')
                r.set_tag('pw', float('%.4g' % v), 'f')
                newq = int(to_mapq(np.array([v]))[0])
                n_changed += int(newq != r.mapping_quality)
                r.mapping_quality = newq
            out.write(r)
    return n_changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True); ap.add_argument('--sd', required=True)
    ap.add_argument('--ref', required=True); ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--batch', type=int, default=50000)
    ap.add_argument('--uncalibrated', action='store_true', help='use the plain GBM instead of GBM + isotonic')
    a = ap.parse_args()
    M = pickle.load(open(a.model, 'rb'))
    if not a.uncalibrated and 'calibrator' in M:
        M['model'] = (M['model_cal'], M['calibrator'])
    sd = SDIndex(read_sd_table(a.sd, a.ref))
    inp = pysam.AlignmentFile(a.inp, 'rb' if a.inp.endswith('.bam') else 'r', check_sq=False)
    hdr = inp.header.to_dict()
    hdr.setdefault('PG', []).append(dict(ID='l5_apply', PN='l5_apply.py', CL=' '.join(sys.argv),
                                         DS='MAPQ recalibrated by %s (sklearn %s)' % (M['run'], M['sklearn'])))
    out = pysam.AlignmentFile(a.out, 'wb', header=hdr)
    batch, cur, recs, n_reads, n_changed = [], None, [], 0, 0
    for r in inp:
        if r.query_name != cur:
            if recs:
                batch.append((cur, recs)); n_reads += 1
                if len(batch) >= a.batch:
                    n_changed += process(batch, M['model'], M['features'], sd, out); batch = []
            cur, recs = r.query_name, []
        recs.append(r)
    if recs:
        batch.append((cur, recs)); n_reads += 1
    if batch:
        n_changed += process(batch, M['model'], M['features'], sd, out)
    out.close()
    sys.stderr.write('reads %d, records with changed MAPQ %d\n' % (n_reads, n_changed))


if __name__ == '__main__':
    main()
