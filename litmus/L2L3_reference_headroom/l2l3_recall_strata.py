#!/usr/bin/env python3
"""L2/L3 pilot: fraction of truth heterozygous SNVs present (same position and ALT base, any
genotype) in each personalization call set, by GRCh38 stratum, for all and 'clean' sites.

  l2l3_recall_strata.py --truth truth.chr20.het_snv.tsv.gz --annot sites/site_annot.tsv \
      --vcf 1kGP_phased=panel.vcf.gz --vcf twopass=twopass.vcf.gz --out recall_by_stratum.tsv
"""
import argparse
import numpy as np
import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--annot', required=True)
    ap.add_argument('--vcf', action='append', required=True)
    ap.add_argument('--chrom', default='chr20')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    het = pd.read_csv(a.truth, sep='\t')
    het['site_id'] = het.chrom + ':' + het.pos0.astype(str)
    het = het.merge(pd.read_csv(a.annot, sep='\t')[['site_id', 'stratum', 'clean']], on='site_id')
    rows = []
    for spec in a.vcf:
        name, path = spec.split('=')
        calls = {}
        for r in pysam.VariantFile(path).fetch(a.chrom):
            if len(r.ref) == 1 and r.alts:
                for alt in r.alts:
                    if len(alt) == 1:
                        calls[r.pos - 1] = alt.upper()
        hit = np.array([calls.get(p) == x for p, x in zip(het.pos0, het.alt)])
        for siteset, m in (('all', np.ones(len(het), bool)), ('clean', het.clean.values == 1)):
            for s in ['CEN', 'SD', 'TR', 'other']:
                mm = m & (het.stratum.values == s)
                rows.append(dict(callset=name, siteset=siteset, stratum=s, truth_het_snv=int(mm.sum()),
                                 present_frac=float(hit[mm].mean()) if mm.any() else np.nan))
    pd.DataFrame(rows).to_csv(a.out, sep='\t', index=False, float_format='%.4f')


if __name__ == '__main__':
    main()
