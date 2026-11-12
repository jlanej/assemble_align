#!/usr/bin/env python3
"""l5_recal.py - MAPQ calibration curves and a learned MAPQ recalibrator (L5 pilot).

For every run (reference x aligner/preset) the evaluable reads (lifted truth valid, read mapped) are
labelled wrong/correct (l5_extract.py). Three predictors of P(wrong) are compared on held-out data:
  RAW  p = 10^(-MAPQ/10) (SAM definition), clipped to [1e-6, 1-1e-6]
  ISO  MAPQ-only recalibration: isotonic regression of Jeffreys-smoothed per-MAPQ error, fit on the
       training fold (what a calibration lookup table can achieve without new information)
  GBM  sklearn HistGradientBoostingClassifier on SAM/PAF-derived + reference-annotation features
  GBMCAL  GBM refit on a random 70% of the training fold + isotonic post-calibration fitted on the other 30%
Chromosome hold-out: truth chromosome parity (odd autosomes + X = fold 1; even, Y, M, unplaced = 0);
train on one fold, predict the other, then swap; metrics use the pooled out-of-fold predictions.
Long reads carry inverse-probability weights from SD-stratified sampling; all statistics are weighted.

usage: l5_recal.py --feat DIR_WITH_*.feat.tsv.gz --outdir OUT [--runs r1,r2,...] [--save-models DIR]
"""
import argparse, glob, json, os, sys, time
for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_v, '2')   # thread budget (set before numpy/sklearn import)
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score, roc_auc_score

EPS = 1e-6
PHRED_BINS = [0, 3, 10, 20, 30, 40, 50, 60, 1000]          # bins on predicted Phred for ECE
MAPQ_BINS = [(0, 0), (1, 3), (4, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 60)]
TARGET_ERR = (1e-3, 1e-4)

BASE_FEATS = ['mapq', 'AS', 'XS', 'ms', 'nn', 'cm', 's1', 's2', 'rl', 'de', 'NM', 'n_sec', 'n_supp', 'n_xa', 'n_sa',
              'clip_l', 'clip_r', 'qlen', 'qalen', 'rspan', 'is_supp', 'proper', 'mate_unmapped', 'tlen',
              'mate_concordant', 'mate_mapq', 'mate_AS', 'mate_XS', 'mate_s2', 'map_sdfrac', 'map_sdid', 'n_nonsec']


def log(*a):
    sys.stderr.write(time.strftime('[%H:%M:%S] ') + ' '.join(map(str, a)) + '\n'); sys.stderr.flush()


def add_derived(d):
    d = d.copy()
    d['as_minus_xs'] = d['AS'] - d['XS']
    d['s2_over_s1'] = d['s2'] / d['s1'].replace(0, np.nan)
    d['sec_as_rel'] = np.where(d['sec_as_max'] >= 0, d['sec_as_max'] / d['AS'].replace(0, np.nan), np.nan)
    d['aln_frac'] = d['qalen'] / d['qlen'].replace(0, np.nan)
    d['nm_rate'] = d['NM'] / d['qalen'].replace(0, np.nan)
    d['rl_frac'] = d['rl'] / d['qlen'].replace(0, np.nan)
    d['clip_frac'] = (d['clip_l'] + d['clip_r']) / d['qlen'].replace(0, np.nan)
    d['span_ratio'] = d['rspan'] / d['qalen'].replace(0, np.nan)
    d['mate_as_minus_xs'] = d.get('mate_AS', np.nan) - d.get('mate_XS', np.nan)
    d['log_tlen'] = np.log10(d['tlen'].clip(lower=1))
    return d


DERIVED = ['as_minus_xs', 's2_over_s1', 'sec_as_rel', 'aln_frac', 'nm_rate', 'rl_frac', 'clip_frac', 'span_ratio',
           'mate_as_minus_xs', 'log_tlen']


def feature_list(d):
    cols = [c for c in BASE_FEATS + DERIVED if c in d.columns and d[c].notna().any() and d[c].nunique(dropna=True) > 1]
    return cols


def p_raw(mapq):
    return np.clip(10 ** (-np.asarray(mapq, float) / 10), EPS, 1 - EPS)


def to_mapq(p):
    return np.clip(np.rint(-10 * np.log10(np.clip(p, EPS, 1))), 0, 60).astype(int)


def wmean(x, w):
    return float(np.sum(x * w) / np.sum(w)) if np.sum(w) > 0 else np.nan


def ece(p, y, w):
    q = -10 * np.log10(np.clip(p, EPS, 1))
    tot = np.sum(w); e = 0.0
    for lo, hi in zip(PHRED_BINS[:-1], PHRED_BINS[1:]):
        m = (q >= lo) & (q < hi)
        if m.any():
            e += np.sum(w[m]) / tot * abs(wmean(p[m], w[m]) - wmean(y[m], w[m]))
    return e


def brier(p, y, w):
    return wmean((p - y) ** 2, w)


def logloss(p, y, w):
    p = np.clip(p, EPS, 1 - EPS)
    return wmean(-(y * np.log(p) + (1 - y) * np.log(1 - p)), w)


def retained_at_error(p, y, w, target):
    """Largest confidence-ordered prefix (ties kept together) whose weighted error <= target.
    returns (weighted correct retained, fraction of weighted correct retained, threshold p)."""
    o = np.argsort(p, kind='stable')
    ps, ys, ws = p[o], y[o], w[o]
    # block boundaries at changes of p
    last = np.r_[np.flatnonzero(np.diff(ps) != 0), len(ps) - 1]
    cw = np.cumsum(ws)[last]; ce = np.cumsum(ws * ys)[last]
    ok = ce / cw <= target
    totc = np.sum(w * (1 - y))
    if not ok.any():
        return 0.0, 0.0, np.nan
    i = np.flatnonzero(ok)[-1]
    corr = cw[i] - ce[i]
    return float(corr), float(corr / totc), float(ps[last[i]])


def precision_at_recall(p, y, w, recall):
    o = np.argsort(p, kind='stable')
    ps, ys, ws = p[o], y[o], w[o]
    last = np.r_[np.flatnonzero(np.diff(ps) != 0), len(ps) - 1]
    cw = np.cumsum(ws)[last]; ce = np.cumsum(ws * ys)[last]
    cc = cw - ce
    totc = np.sum(w * (1 - y))
    ok = cc / totc >= recall
    if not ok.any():
        return np.nan, np.nan
    i = np.flatnonzero(ok)[0]
    return float(cc[i] / cw[i]), float(cc[i] / totc)


def selfselect(p, y, w, target):
    """Label-free selection: largest confidence-ordered prefix (ties together) whose mean *predicted*
    error is <= target; returns realised weighted error and fraction of correct reads retained."""
    o = np.argsort(p, kind='stable')
    ps, ys, ws = p[o], y[o], w[o]
    last = np.r_[np.flatnonzero(np.diff(ps) != 0), len(ps) - 1]
    cw = np.cumsum(ws)[last]; cp = np.cumsum(ws * ps)[last]; ce = np.cumsum(ws * ys)[last]
    ok = cp / cw <= target
    totc = np.sum(w * (1 - y))
    if not ok.any():
        return np.nan, 0.0
    i = np.flatnonzero(ok)[-1]
    return float(ce[i] / cw[i]), float((cw[i] - ce[i]) / totc)


def fit_iso(mapq, y, w):
    df = pd.DataFrame(dict(q=mapq, y=y * w, w=w)).groupby('q').sum()
    pj = (df.y + 0.5) / (df.w + 1.0)             # Jeffreys-smoothed per-MAPQ error
    ir = IsotonicRegression(increasing=False, y_min=EPS, y_max=1 - EPS, out_of_bounds='clip')
    ir.fit(df.index.values.astype(float), pj.values, sample_weight=df.w.values)
    return ir


def fit_gbm(X, y, w, seed=0):
    mono = [(-1 if c == 'mapq' else 0) for c in X.columns]
    m = HistGradientBoostingClassifier(learning_rate=0.05, max_iter=600, max_leaf_nodes=31,
                                       min_samples_leaf=40, l2_regularization=1.0, early_stopping=True,
                                       validation_fraction=0.15, n_iter_no_change=40, monotonic_cst=mono,
                                       random_state=seed)
    m.fit(X, y, sample_weight=w)
    return m


def run_one(path, outdir, save_dir=None):
    d = pd.read_csv(path, sep='\t', low_memory=False)
    run = d.run.iloc[0]
    n_all = len(d); n_unmapped = int((d.mapped == 0).sum()); n_notruth = int((d.truth_ok == 0).sum())
    d = d[(d.mapped == 1) & (d.truth_ok == 1) & (d.fold >= 0)].copy()
    d = add_derived(d)
    feats = feature_list(d)
    y = (1 - d.correct.values).astype(int)
    w = d.weight.values.astype(float)
    mq = d.mapq.values.astype(float)
    pr = p_raw(mq)
    piso = np.full(len(d), np.nan); pgbm = np.full(len(d), np.nan); pgc = np.full(len(d), np.nan)
    X = d[feats]
    imps = []
    for f in (0, 1):
        tr, te = d.fold.values != f, d.fold.values == f
        ir = fit_iso(mq[tr], y[tr], w[tr])
        piso[te] = ir.predict(mq[te])
        if y[tr].sum() >= 5:
            m = fit_gbm(d.loc[tr, feats], y[tr], w[tr])
            pgbm[te] = np.clip(m.predict_proba(d.loc[te, feats])[:, 1], EPS, 1 - EPS)
            imps.append(dict(fold_test=f, n_iter=int(m.n_iter_)))
        else:
            pgbm[te] = piso[te]
        # GBM + held-out isotonic post-calibration (fit GBM on 70% of the training fold, isotonic on 30%)
        idx = np.flatnonzero(tr)
        sel = np.random.default_rng(11 + f).random(len(idx)) < 0.7
        A, B = idx[sel], idx[~sel]
        if y[A].sum() >= 5 and y[B].sum() >= 3:
            mA = fit_gbm(X.iloc[A], y[A], w[A])
            irc = IsotonicRegression(increasing=True, y_min=EPS, y_max=1 - EPS, out_of_bounds='clip')
            irc.fit(mA.predict_proba(X.iloc[B])[:, 1], y[B], sample_weight=w[B])
            pgc[te] = np.clip(irc.predict(mA.predict_proba(X.loc[te])[:, 1]), EPS, 1 - EPS)
        else:
            pgc[te] = pgbm[te]
    if save_dir and y.sum() >= 5:
        import pickle
        os.makedirs(save_dir, exist_ok=True)
        mfull = fit_gbm(d[feats], y, w)
        sel = np.random.default_rng(7).random(len(d)) < 0.7
        mcal = fit_gbm(d[feats].iloc[np.flatnonzero(sel)], y[sel], w[sel])
        cal = IsotonicRegression(increasing=True, y_min=EPS, y_max=1 - EPS, out_of_bounds='clip')
        cal.fit(mcal.predict_proba(d[feats].iloc[np.flatnonzero(~sel)])[:, 1], y[~sel], sample_weight=w[~sel])
        with open(os.path.join(save_dir, 'model_%s.pkl' % run), 'wb') as fh:
            pickle.dump(dict(model=mfull, model_cal=mcal, calibrator=cal, features=feats, run=run,
                             sklearn=__import__('sklearn').__version__), fh)
    d['p_raw'], d['p_iso'], d['p_gbm'], d['p_gbmcal'] = pr, piso, pgbm, pgc
    d['mapq_recal'] = to_mapq(pgbm)
    d['y_wrong'] = y
    keep = ['qname', 'mate', 'run', 'fold', 'in_sd', 'ambig', 'weight', 'mapq', 'mapq_recal', 'p_raw', 'p_iso', 'p_gbm', 'p_gbmcal',
            'y_wrong', 'map_sdfrac', 't_sdid', 't_sdfrac']
    d[keep].to_csv(os.path.join(outdir, 'pred_%s.tsv.gz' % run), sep='\t', index=False, compression='gzip')
    meta = dict(run=run, n_reads=n_all, n_unmapped=n_unmapped, n_no_truth=n_notruth, n_eval=len(d),
                features=feats, gbm=imps)
    return d, meta


def metrics_rows(d, run):
    rows = []
    for stratum, m in (('all', np.ones(len(d), bool)), ('in_SD', d.in_sd.values == 1), ('outside_SD', d.in_sd.values == 0),
                       ('all_excl_ambig_truth', d.ambig.values == 0)):
        if m.sum() == 0:
            continue
        y, w = d.y_wrong.values[m], d.weight.values[m]
        for pred in [x for x in ('raw', 'iso', 'gbm', 'gbmcal') if 'p_' + x in d.columns]:
            p = d['p_' + pred].values[m]
            r = dict(run=run, stratum=stratum, predictor=pred, n_reads=int(m.sum()), n_wrong=int(y.sum()),
                     w_error_rate=wmean(y, w), ECE=ece(p, y, w), Brier=brier(p, y, w), logloss=logloss(p, y, w),
                     AUPRC_wrong=(average_precision_score(y, p, sample_weight=w) if 0 < y.sum() < len(y) else np.nan),
                     AUROC=(roc_auc_score(y, p, sample_weight=w) if 0 < y.sum() < len(y) else np.nan))
            for t in TARGET_ERR:
                c, fr, thr = retained_at_error(p, y, w, t)
                r['correct_retained_err%.0e' % t] = c; r['frac_correct_retained_err%.0e' % t] = fr
            for t in TARGET_ERR:
                re_, fr_ = selfselect(p, y, w, t)
                r['selfsel_realised_err%.0e' % t] = re_; r['selfsel_frac_correct%.0e' % t] = fr_
            qq = d.mapq.values[m] if pred == 'raw' else to_mapq(p)
            keep30 = qq >= 30
            r['mapq30_realised_err'] = wmean(y[keep30], w[keep30]) if keep30.any() else np.nan
            r['mapq30_frac_correct'] = float(np.sum(w[keep30] * (1 - y[keep30])) / np.sum(w * (1 - y)))
            for rc in (0.99, 0.995):
                pr_, rr = precision_at_recall(p, y, w, rc)
                r['precision_at_recall%.3f' % rc] = pr_
            q = d.mapq.values[m]
            r['n_mapq60'] = int((q == 60).sum()); r['n_wrong_mapq60'] = int(((q == 60) & (y == 1)).sum())
            r['w_err_mapq60'] = wmean(y[q == 60], w[q == 60]) if (q == 60).any() else np.nan
            r['n_mapq0'] = int((q == 0).sum()); r['n_wrong_mapq0'] = int(((q == 0) & (y == 1)).sum())
            r['w_err_mapq0'] = wmean(y[q == 0], w[q == 0]) if (q == 0).any() else np.nan
            rows.append(r)
    return rows


def jeffreys_ci(k, n):
    from scipy.stats import beta
    lo = beta.ppf(0.025, k + 0.5, n - k + 0.5) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 0.5, n - k + 0.5) if k < n else 1.0
    return lo, hi


def calib_rows(d, run):
    rows = []
    for binned_by, qv in (('reported_mapq', d.mapq.values), ('recalibrated_mapq', d.mapq_recal.values)):
      for stratum, m in (('all', np.ones(len(d), bool)), ('in_SD', d.in_sd.values == 1), ('outside_SD', d.in_sd.values == 0)):
        for lo, hi in MAPQ_BINS:
            b = m & (qv >= lo) & (qv <= hi)
            n = int(b.sum())
            if n == 0:
                continue
            y, w = d.y_wrong.values[b], d.weight.values[b]
            k = int(y.sum())
            ew = wmean(y, w)
            neff = (w.sum() ** 2) / (w ** 2).sum()          # Kish effective sample size
            lo_, hi_ = jeffreys_ci(ew * neff, neff)
            rows.append(dict(run=run, binned_by=binned_by, stratum=stratum, mapq_bin='%d-%d' % (lo, hi) if lo != hi else str(lo), n_reads=n,
                             n_wrong=k, w_error=ew, err_ci_lo=lo_, err_ci_hi=hi_,
                             empirical_phred=(-10 * np.log10(ew) if ew > 0 else np.nan),
                             empirical_phred_lower_bound=(-10 * np.log10(hi_) if hi_ > 0 else np.nan),
                             mean_reported_p=wmean(d.p_raw.values[b], w), mean_gbm_p=wmean(d.p_gbm.values[b], w),
                             mean_iso_p=wmean(d.p_iso.values[b], w),
                             mean_gbmcal_p=(wmean(d.p_gbmcal.values[b], w) if 'p_gbmcal' in d.columns else np.nan)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--feat', required=True); ap.add_argument('--outdir', required=True)
    ap.add_argument('--runs', default='')
    ap.add_argument('--save-models', default='', help='dir for final models trained on all evaluable reads')
    ap.add_argument('--metrics-only', action='store_true', help='recompute metrics from OUTDIR/pred_*.tsv.gz')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.feat, '*.feat.tsv.gz')))
    if a.runs:
        keep = set(a.runs.split(','))
        files = [f for f in files if os.path.basename(f).split('.feat')[0] in keep]
    met, cal, metas = [], [], []
    if a.metrics_only:
        for f in sorted(glob.glob(os.path.join(a.outdir, 'pred_*.tsv.gz'))):
            d = pd.read_csv(f, sep='\t', low_memory=False)
            run = d.run.iloc[0]
            met += metrics_rows(d, run); cal += calib_rows(d, run)
            log('metrics', run)
        files = []
    for f in files:
        log('run', f)
        d, meta = run_one(f, a.outdir, a.save_models or None)
        met += metrics_rows(d, meta['run']); cal += calib_rows(d, meta['run']); metas.append(meta)
        log(meta['run'], 'eval', meta['n_eval'], 'wrong', int(d.y_wrong.sum()), 'gbm', meta['gbm'])
    pd.DataFrame(met).to_csv(os.path.join(a.outdir, 'L5_recalibrator_metrics.tsv'), sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(cal).to_csv(os.path.join(a.outdir, 'L5_calibration.tsv'), sep='\t', index=False, float_format='%.6g')
    if metas:
        json.dump(metas, open(os.path.join(a.outdir, 'L5_run_meta.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
