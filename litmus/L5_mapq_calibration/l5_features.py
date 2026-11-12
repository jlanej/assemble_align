#!/usr/bin/env python3
"""l5_features.py - per-read alignment features shared by l5_extract.py (training/evaluation) and
l5_apply.py (rewriting MAPQ). Record selection follows paftools.js mapeval mode 0: among
non-secondary mapped records the one with most matching bases is the evaluated record."""
import numpy as np
import pandas as pd

TAGS_I = ['NM', 'AS', 'XS', 'ms', 'nn', 'cm', 's1', 's2', 'rl']
TAGS_F = ['de']


def rec_stats(r):
    mlen = gap = clip_l = clip_r = 0
    seen = False
    for op, ln in (r.cigartuples or []):
        if op in (0, 7, 8):
            mlen += ln; seen = True
        elif op in (1, 2):
            gap += ln; seen = True
        elif op in (4, 5):
            if seen:
                clip_r += ln
            else:
                clip_l += ln
    nm = r.get_tag('NM') if r.has_tag('NM') else None
    score = (mlen - (nm - gap)) if (nm is not None and nm >= gap) else (r.reference_end - r.reference_start)
    return score, mlen, gap, clip_l, clip_r


def new_state():
    return dict(best=None, best_score=-1, n_sec=0, n_supp=0, sec_as_max=-1, qlen=0, n_nonsec=0, best_id=None)


def mate_of(r):
    return 1 if r.is_read1 else (2 if r.is_read2 else 0)


def update(d, r, rec_id=None):
    if r.query_sequence is not None and not r.is_secondary and not r.is_supplementary:
        d['qlen'] = r.infer_read_length() or len(r.query_sequence)
    if r.is_unmapped:
        return
    if r.is_secondary:
        d['n_sec'] += 1
        if r.has_tag('AS'):
            d['sec_as_max'] = max(d['sec_as_max'], r.get_tag('AS'))
        return
    if r.is_supplementary:
        d['n_supp'] += 1
    d['n_nonsec'] += 1
    score, mlen, gap, cl, cr = rec_stats(r)
    if score > d['best_score']:
        f = dict(chrom=r.reference_name, pos=r.reference_start, end=r.reference_end,
                 strand=-1 if r.is_reverse else 1, mapq=min(60, r.mapping_quality),
                 is_supp=int(r.is_supplementary), mlen=mlen, gap=gap, clip_l=cl, clip_r=cr,
                 qalen=r.query_alignment_length, rspan=r.reference_end - r.reference_start,
                 proper=int(r.is_proper_pair), mate_unmapped=int(r.mate_is_unmapped) if r.is_paired else 0,
                 tlen=abs(r.template_length),
                 mate_same_chr=int(r.is_paired and not r.mate_is_unmapped and r.next_reference_id == r.reference_id),
                 n_xa=len(r.get_tag('XA').rstrip(';').split(';')) if r.has_tag('XA') else 0,
                 n_sa=len(r.get_tag('SA').rstrip(';').split(';')) if r.has_tag('SA') else 0,
                 tp=r.get_tag('tp') if r.has_tag('tp') else '')
        for t in TAGS_I:
            f[t] = r.get_tag(t) if r.has_tag(t) else np.nan
        for t in TAGS_F:
            f[t] = r.get_tag(t) if r.has_tag(t) else np.nan
        d['best'] = f; d['best_score'] = score; d['best_id'] = rec_id


def rows_from(states):
    rows = []
    for (q, m), d in states.items():
        row = dict(qname=q, mate=m, qlen=d['qlen'], n_sec=d['n_sec'], n_supp=d['n_supp'],
                   sec_as_max=d['sec_as_max'], n_nonsec=d['n_nonsec'], mapped=int(d['best'] is not None))
        if d['best'] is not None:
            row.update(d['best'])
        rows.append(row)
    return pd.DataFrame(rows)


def add_mate_and_sd(df, sd):
    if (df.mate > 0).any():
        cols = ['qname', 'mate', 'mapq', 'AS', 'XS', 's2', 'chrom', 'pos']
        for c in cols:
            if c not in df.columns:
                df[c] = np.nan
        mt = df[cols].copy()
        mt['mate'] = 3 - mt['mate']
        mt.columns = ['qname', 'mate', 'mate_mapq', 'mate_AS', 'mate_XS', 'mate_s2', 'mate_chrom', 'mate_pos']
        df = df.merge(mt, on=['qname', 'mate'], how='left')
        df['mate_concordant'] = ((df.mate_chrom == df.chrom) & ((df.mate_pos - df.pos).abs() < 2000)).astype(int)
        df = df.drop(columns=['mate_chrom', 'mate_pos'])
    df['map_sdfrac'] = 0.0; df['map_sdid'] = 0.0
    mp = df.mapped.values == 1
    if mp.any():
        f_, i_ = sd.query(df.chrom.values[mp], df.pos.values[mp].astype(np.int64), df.end.values[mp].astype(np.int64))
        df.loc[mp, 'map_sdfrac'] = f_; df.loc[mp, 'map_sdid'] = i_
    return df
