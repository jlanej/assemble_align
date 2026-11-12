#!/usr/bin/env python3
"""L2/L3 pilot: render L2L3_reference_headroom.md from the result tables (all numbers are read
from the TSVs so that the report and the tables cannot drift apart).

  l2l3_report.py --results results/ --links links.json --out L2L3_reference_headroom.md
links.json maps artifact file names to markdown link targets (optional).
"""
import argparse, json, os
import numpy as np
import pandas as pd

REFS = ['GRCh38', 'CHM13', 'pers2p', 'pers1KG', 'oracle']
RNAME = {'GRCh38': 'GRCh38', 'CHM13': 'CHM13', 'pers2p': 'GRCh38+own calls (unphased)',
         'pers1KG': 'GRCh38+1kGP phased', 'oracle': 'own assembly (oracle)'}
TNAME = {'ilmn': 'Illumina', 'hifi': 'HiFi', 'ont': 'ONT'}


def pct(x, d=2):
    return ('%.' + str(d) + 'f%%') % (100 * x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)
    ap.add_argument('--links', default=None)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tail', default=None, help='static markdown appended after the generated sections')
    a = ap.parse_args()
    R = a.results
    J = json.load(open(a.links)) if a.links else {}
    L = J.get('links', J) if isinstance(J, dict) else {}
    D = J.get('desc', {}) if isinstance(J, dict) else {}
    lk = lambda f: '[%s](%s)' % (f, L.get(f, f))
    m = pd.read_csv(os.path.join(R, 'L2L3_metrics.tsv'), sep='\t').set_index(['tech', 'ref'])
    st = pd.read_csv(os.path.join(R, 'L2L3_metrics_by_stratum.tsv'), sep='\t').set_index(['tech', 'stratum', 'ref'])
    ci = pd.read_csv(os.path.join(R, 'L2L3_readmetric_ci.tsv'), sep='\t')
    abci = pd.read_csv(os.path.join(R, 'L2L3_allele_balance_ci.tsv'), sep='\t')
    ab = pd.read_csv(os.path.join(R, 'L2L3_allele_balance.tsv'), sep='\t')
    ab = ab[ab.mapq_min.astype(str) == '20']
    cs = pd.read_csv(os.path.join(R, 'L2L3_callsets_vs_truth.tsv'), sep='\t').set_index('callset')
    tr = pd.read_csv(os.path.join(R, 'L2L3_mapq_transitions.tsv'), sep='\t').set_index(['tech', 'ref'])
    truth = json.load(open(os.path.join(R, 'L2L3_truth_summary.json')))
    ctx = {t: pd.read_csv(os.path.join(R, 'L2L3_lr_error_context_%s.tsv' % t), sep='\t').set_index('context') for t in ('hifi', 'ont')}
    qv = {t: pd.read_csv(os.path.join(R, 'L2L3_lr_qv_calibration_%s.tsv' % t), sep='\t') for t in ('hifi', 'ont')}

    def C(tech, stratum, ref, metric='frac_mapq_lt20'):
        x = ci[(ci.tech == tech) & (ci.stratum == stratum) & (ci.ref == ref) & (ci.metric == metric)].iloc[0]
        return x

    def rec(tech, stratum, ref, metric='frac_mapq_lt20'):
        x = C(tech, stratum, ref, metric)
        return '%.0f%% [%.0f, %.0f]' % (100 * x.gap_recovered, 100 * x.gap_recovered_ci_lo, 100 * x.gap_recovered_ci_hi)

    def AB(tech, stratum, ref, siteset='clean'):
        return abci[(abci.tech == tech) & (abci.stratum == stratum) & (abci.ref == ref) & (abci.siteset == siteset)].iloc[0]

    def abtxt(tech, stratum, ref):
        x = AB(tech, stratum, ref)
        return '%.3f [%.3f, %.3f]' % (x.ref_frac, x.ci_lo, x.ci_hi)

    def abrec(tech, stratum, ref):
        g, o, x = AB(tech, stratum, 'GRCh38').ref_frac, AB(tech, stratum, 'oracle').ref_frac, AB(tech, stratum, ref).ref_frac
        return 100 * (g - x) / (g - o) if g != o else float('nan')

    def drop(tech, stratum, ref):
        x = ab[(ab.tech == tech) & (ab.stratum == stratum) & (ab.ref == ref) & (ab.siteset == 'clean')].iloc[0]
        return x.frac_sites_alt0_depth10

    def gapfrac(tech, stratum, col, ref, src=st):
        if stratum is None:
            g, o, v = src.loc[(tech, 'GRCh38'), col], src.loc[(tech, 'oracle'), col], src.loc[(tech, ref), col]
        else:
            g, o, v = src.loc[(tech, stratum, 'GRCh38'), col], src.loc[(tech, stratum, 'oracle'), col], src.loc[(tech, stratum, ref), col]
        return 100 * (g - v) / (g - o)

    n = {t: int(m.loc[(t, 'GRCh38'), 'n_reads']) for t in ('ilmn', 'hifi', 'ont')}
    sr = {s: int(st.loc[('ilmn', s, 'GRCh38'), 'n_reads']) for s in ('CEN', 'SD', 'TR', 'other', 'unplaced')}
    lrs = {t: {s: int(st.loc[(t, s, 'GRCh38'), 'n_reads']) for s in ('CEN', 'SD', 'nonSD', 'unplaced')} for t in ('hifi', 'ont')}
    nsite = {s: int(AB('ilmn', s, 'GRCh38').n_sites) for s in ('ALL', 'SD', 'TR', 'other', 'CEN')}
    nsite_all = {s: int(AB('ilmn', s, 'GRCh38', 'all').n_sites) for s in ('ALL', 'SD', 'TR', 'other', 'CEN')}
    out = []
    w = out.append
    w('# L2/L3 litmus pilot: reference headroom on real 1kG data (NA19240, chr20)\n')
    w('*Track L2/L3, run 2026-11-11. One sample (NA19240, YRI), one chromosome (chr20), long reads at ~3x; '
      'every number below is provisional pilot evidence, not a population estimate.*\n')
    w('## Bottom line\n')
    g, o = m.loc[('ilmn', 'GRCh38'), 'frac_mapq_lt20'], m.loc[('ilmn', 'oracle'), 'frac_mapq_lt20']
    g2, o2 = C('ilmn', 'all_noCEN', 'GRCh38').value, C('ilmn', 'all_noCEN', 'oracle').value
    w('1. **Short reads, mapping confidence.** %s of the %s Illumina reads in the NYGC chr20 slice are unmapped or MAPQ<20 on GRCh38 '
      '(minimap2 `-x sr`); mapping the same reads to the sample\'s own diploid assembly (per-template best haplotype) lowers this to %s. '
      'A realistic two-pass personalized reference (bcftools calls from the same reads, unphased diploid pair) recovers %s of that gap '
      'and a 1kGP-phased personalized diploid reference %s (95%% block-bootstrap CI); CHM13 recovers %s. Outside the chr20 centromere model '
      '(+-1 Mbp) the gap is %s -> %s and the recoveries are %s and %s.' % (
          pct(g), '{:,}'.format(n['ilmn']), pct(o), rec('ilmn', 'all', 'pers2p'), rec('ilmn', 'all', 'pers1KG'), rec('ilmn', 'all', 'CHM13'),
          pct(g2), pct(o2), rec('ilmn', 'all_noCEN', 'pers2p'), rec('ilmn', 'all_noCEN', 'pers1KG')))
    w('2. **Short reads, reference bias.** At %s clean truth heterozygous SNVs (assembly-derived, outside the centromere model, no indel '
      'within 5-10 bp) the GRCh38-REF allele fraction among MAPQ>=20 Illumina reads is %s on GRCh38 and %s on the oracle. Small-variant '
      'personalization removes %.0f%% (own calls) and %.0f%% (1kGP phased) of the excess; CHM13 removes %.0f%%. In segmental duplications '
      '(%s sites) bias falls from %s to %s with the oracle, but personalization removes only %.0f%% / %.0f%% of it, and the fraction of '
      'het SNVs with zero ALT reads at depth>=10 is %s (own calls) and %s (1kGP phased) versus %s on GRCh38 and %s on the oracle.' % (
          '{:,}'.format(nsite['ALL']), abtxt('ilmn', 'ALL', 'GRCh38'), abtxt('ilmn', 'ALL', 'oracle'), abrec('ilmn', 'ALL', 'pers2p'),
          abrec('ilmn', 'ALL', 'pers1KG'), abrec('ilmn', 'ALL', 'CHM13'), '{:,}'.format(nsite['SD']), abtxt('ilmn', 'SD', 'GRCh38'),
          abtxt('ilmn', 'SD', 'oracle'), abrec('ilmn', 'SD', 'pers2p'), abrec('ilmn', 'SD', 'pers1KG'), pct(drop('ilmn', 'SD', 'pers2p'), 1),
          pct(drop('ilmn', 'SD', 'pers1KG'), 1), pct(drop('ilmn', 'SD', 'GRCh38'), 1), pct(drop('ilmn', 'SD', 'oracle'), 1)))
    w('3. **Where the short-read headroom is.** In unique sequence ("other" stratum, %s reads) personalization captures most of the '
      'small-variant headroom: perfect alignments (NM=0, unclipped) rise from %s (GRCh38) to %s (1kGP phased) versus %s (oracle), i.e. '
      '%.0f%% of the gap. Soft-clipping (>20 bp) is structural: CHM13 recovers %s of the clipping gap across all reads, small-variant '
      'personalization %s.' % (
          '{:,}'.format(sr['other']), pct(st.loc[('ilmn', 'other', 'GRCh38'), 'frac_perfect'], 1), pct(st.loc[('ilmn', 'other', 'pers1KG'), 'frac_perfect'], 1),
          pct(st.loc[('ilmn', 'other', 'oracle'), 'frac_perfect'], 1), gapfrac('ilmn', 'other', 'frac_perfect', 'pers1KG'),
          rec('ilmn', 'all', 'CHM13', 'frac_clipped'), rec('ilmn', 'all', 'pers1KG', 'frac_clipped')))
    for t in ('hifi', 'ont'):
        w('%d. **%s (%s chr20 reads, ~3x).** MAPQ<20/unmapped falls from %s (GRCh38) to %s (oracle); personalization recovers %s (own '
          'calls) and %s (1kGP phased) of the gap, CHM13 %s. Primary alignments clipped >1 kb: %s -> %s; CHM13 recovers %s, personalization '
          '%s. Personalization does remove most of the per-read divergence: median gap-compressed divergence %s -> %s (oracle), %s with '
          '1kGP phased (%.0f%% of the gap).' % (
              4 if t == 'hifi' else 5, TNAME[t], '{:,}'.format(n[t]), pct(m.loc[(t, 'GRCh38'), 'frac_mapq_lt20']), pct(m.loc[(t, 'oracle'), 'frac_mapq_lt20']),
              rec(t, 'all', 'pers2p'), rec(t, 'all', 'pers1KG'), rec(t, 'all', 'CHM13'), pct(m.loc[(t, 'GRCh38'), 'frac_clipped']),
              pct(m.loc[(t, 'oracle'), 'frac_clipped']), rec(t, 'all', 'CHM13', 'frac_clipped'), rec(t, 'all', 'pers1KG', 'frac_clipped'),
              pct(m.loc[(t, 'GRCh38'), 'gc_divergence_median'], 2), pct(m.loc[(t, 'oracle'), 'gc_divergence_median'], 2),
              pct(m.loc[(t, 'pers1KG'), 'gc_divergence_median'], 2), gapfrac(t, None, 'gc_divergence_median', 'pers1KG', m)))
    w('6. **Long-read reference bias is concentrated in tandem repeats (ONT).** GRCh38-REF fraction at clean het SNVs in TRs: ONT %s '
      '(oracle %s; own calls %s, 1kGP phased %s, CHM13 %s); HiFi %s (oracle %s). Genome-wide (non-CEN) ONT shows %s versus %s on the oracle.' % (
          abtxt('ont', 'TR', 'GRCh38'), abtxt('ont', 'TR', 'oracle'), abtxt('ont', 'TR', 'pers2p'), abtxt('ont', 'TR', 'pers1KG'),
          abtxt('ont', 'TR', 'CHM13'), abtxt('hifi', 'TR', 'GRCh38'), abtxt('hifi', 'TR', 'oracle'), abtxt('ont', 'ALL', 'GRCh38'), abtxt('ont', 'ALL', 'oracle')))
    hc, oc = ctx['hifi'], ctx['ont']
    err = lambda d: d.mismatch_bases + d.ins_bases + d.del_bases
    hp4 = ['HP4-5', 'HP6-7', 'HP8-10', 'HP11-15', 'HP16+']

    def qbin(t, lo, hi):
        d = qv[t][(qv[t].reported_qv >= lo) & (qv[t].reported_qv < hi)]
        return -10 * np.log10(max(d.errors_mm_or_ins.sum(), 0.5) / d.bases.sum())
    w('7. **Residual long-read errors versus the own assembly.** HiFi (DeepConsensus) %.2f errors/kb overall, rising from %.2f/kb outside '
      'homopolymers to %.1f/kb in runs >=16 bp; runs >=4 bp hold %.1f%% of bases but %.0f%% of error bases, and %.0f%% of error bases are '
      'indels. ONT R10.4.1 (dorado sup 4.1.0, unfiltered) %.1f errors/kb, %.1f/kb outside homopolymers to %.0f/kb in runs >=16 bp. Reported '
      'base QVs are over-confident at the top: HiFi bases reported Q40-49 are empirically Q%.1f and ONT Q40-49 bases Q%.1f '
      '(mismatch+insertion errors; the own-assembly floor limits empirical QV above ~Q40).\n' % (
          hc.loc['ALL', 'err_per_kb'], hc.loc['HP1', 'err_per_kb'], hc.loc['HP16+', 'err_per_kb'], 100 * hc.loc[hp4, 'frac_of_bases'].sum(),
          100 * err(hc)[hp4].sum() / err(hc)['ALL'], 100 * (hc.loc['ALL', 'ins_bases'] + hc.loc['ALL', 'del_bases']) / err(hc)['ALL'],
          oc.loc['ALL', 'err_per_kb'], oc.loc['HP1', 'err_per_kb'], oc.loc['HP16+', 'err_per_kb'], qbin('hifi', 40, 50), qbin('ont', 40, 50)))

    w('![Figure 1](%s)\n' % L.get('L2L3_fig1_reference_headroom.png', 'L2L3_fig1_reference_headroom.png'))
    w('*Figure 1. Reference headroom for NA19240 chr20 reads by technology (columns). Rows: reads unmapped or MAPQ<20; pooled edit '
      'rate excluding indels >=50 bp; GRCh38-REF allele fraction at clean truth heterozygous SNVs (MAPQ>=20) by GRCh38 stratum '
      '(dashed line = unbiased 0.5). Bars: GRCh38, CHM13, GRCh38 personalized with the sample\'s own unphased bcftools calls, GRCh38 '
      'personalized with the 1kGP phased panel, and the sample\'s own HPRC release-2 assembly (oracle). Rows 1-2 include centromere-model '
      'reads (13.5% of Illumina reads), which dominate the totals; strata are in the tables.*\n')

    w('## Design\n')
    w('The pilot asks two questions on public 1kG-linked data. **L2 (oracle headroom):** how much better do the same reads align when '
      'the reference *is* the sample, i.e. the HPRC release-2 (preprint; [Lucas 2026](https://doi.org/10.64898/2026.07.21.739710)) trio-phased hifiasm '
      '([Cheng 2021](https://doi.org/10.1038/s41592-020-01056-5)) assembly of NA19240 (hap1 = paternal '
      '`NA19240_pat_hprc_r2_v1.0.1`, hap2 = maternal), mapping every read to each whole haplotype separately and keeping, per template, '
      'the haplotype with the larger summed alignment score? **L3 (realistic recovery):** what fraction of that gap does a method that '
      'could be run on any 1kGP sample recover? Two realistic arms were built on chr20 of GRCh38 (the rest of the genome left as GRCh38): '
      '(i) **two-pass**, bcftools mpileup/call ([Danecek 2021](https://doi.org/10.1093/gigascience/giab008)) on the GRCh38 minimap2 alignment of the same Illumina reads (QUAL>=30, 8<=DP<=90), applied as '
      'an unphased diploid pair (A = all ALT alleles, R = homozygous ALT only); (ii) **1kGP phased**, the sample\'s genotypes from the 1kGP '
      'high-coverage phased SNV/INDEL/SV panel ([Byrska-Bishop 2022](https://doi.org/10.1016/j.cell.2022.08.004); sequence-resolved alleles; 137 symbolic SV records dropped) applied as haplotype 1 and 2. '
      'All references are whole genomes (no single-chromosome indexes); the same minimap2 2.31 ([Li 2018](https://doi.org/10.1093/bioinformatics/bty191)) preset is used for every reference: '
      '`-x sr` (Illumina), `-x map-hifi` (HiFi), `-x lr:hq` (ONT R10.4.1), with `--eqx -Y --secondary=no` for long reads.\n')
    w('| Read set | Source | Selection | Reads | Notes |\n|---|---|---|---|---|')
    w('| Illumina 2x150, 30x | NYGC CRAM ERR3989410 (1kGP 698-related collection; md5 verified) | all primary, non-duplicate pairs with >=1 record on chr20 or mate on chr20 (mates elsewhere and unmapped mates included) | %s reads (%s pairs) | selection by NYGC bwa-mem GRCh38 placement |' % ('{:,}'.format(n['ilmn']), '{:,}'.format(n['ilmn'] // 2)))
    w('| PacBio HiFi | HPRC r2 DeepConsensus v1.2 q20 FASTQ, 3 SMRT cells | first 3.1 Gbp of each file (9.3 Gbp, ~3x), then chr20 prefilter + union rule | %s reads | DC reads were also used to build the assembly |' % '{:,}'.format(n['hifi']))
    w('| ONT R10.4.1 UL | HPRC r2 dorado 0.7.2 sup@v4.1.0 BAMs, 3 runs (read N50 ~100 kb) | first 3.1 Gbp of each file (9.3 Gbp, ~3x), then chr20 prefilter + union rule | %s reads | no read-QV filter applied |' % '{:,}'.format(n['ont']))
    w('')
    w('**Long-read selection.** Each streamed read was first mapped without base-level alignment to an index of the four chr20 sequences '
      '(GRCh38, CHM13, hap1 `CM099588.1`, hap2 `CM099626.1`); reads whose best single chain covered >=80% of the read were candidates '
      '(HiFi 17,726; ONT 6,184). Candidates were mapped to every whole-genome reference, and the analysed set is the union of reads whose '
      'primary alignment lies on chr20 in GRCh38, CHM13 or the oracle, so no single reference defines the set.\n')
    w('**Truth.** hap1/hap2 chr20 were aligned to whole-genome GRCh38 and CHM13 (`minimap2 -cx asm5 --cs -z200000,10000 -r2k`); '
      'dipcall-like rules ([Li 2018](https://doi.org/10.1038/s41592-018-0054-7)) kept GRCh38 chr20 bases covered by exactly one >=10 kb alignment of each haplotype (%s bp confident of %s ACGT bp), '
      'giving %s heterozygous and %s homozygous-ALT SNVs; each het SNV carries its coordinate on hap1, hap2, CHM13 (both haplotype projections '
      'must agree; %s sites) and on every personalized haplotype (bcftools consensus chain). The GRCh38 centromere model +-1 Mbp '
      '(chr20:25.44-31.04 Mb; %s het sites, %s of them clean) is reported separately because haplotype alpha-satellite aligned to the model sequence is not '
      'orthology. "Clean" sites additionally exclude SNVs within 5 bp of a truth indel or within 10 bp of any indel/complex record of either '
      'personalization call set, because equivalent haplotypes written with different indel/SNV combinations otherwise produce artefactual '
      'allele counts at a single coordinate.\n' % (
          '{:,}'.format(truth['confident_bp']), '{:,}'.format(truth['acgt_bp']), '{:,}'.format(truth['het_snv']), '{:,}'.format(truth['homalt_snv']),
          '{:,}'.format(truth['het_snv_with_chm13']), '{:,}'.format(nsite_all['CEN']), '{:,}'.format(nsite['CEN'])))
    w('**Strata** come from each read\'s GRCh38 primary alignment so that every reference is scored on identical read sets: CEN (>=50% of '
      'the aligned span in a centromere model +-1 Mbp), SD (>=50% in UCSC genomicSuperDups), TR (any overlap with UCSC simpleRepeat or '
      'homopolymers >=7 bp, padded 5 bp; short reads only), other, and unplaced (unmapped on GRCh38). Sites use the same intervals by '
      'position. **Statistics:** 95% intervals are block bootstraps (1,000 replicates, 1 Mbp GRCh38 blocks, identical weights for all '
      'references within a technology, so differences and gap recoveries are paired). Gap recovered = (GRCh38 - X)/(GRCh38 - oracle).\n')

    w('## Results\n')
    w('### Read-level metrics, all reads\n')
    cols = [('frac_unmapped', 'unmapped', 2), ('frac_mapq0', 'MAPQ0', 2), ('frac_mapq_lt20', 'MAPQ<20 or unmapped', 2),
            ('frac_clipped', 'clipped (>20 bp SR, >1 kb LR)', 2), ('frac_with_supp', 'with supplementary', 2),
            ('edit_rate_small_pooled', 'edit rate (indels<50 bp)', 3)]
    w('| Tech | Reference | ' + ' | '.join(c[1] for c in cols) + ' | proper pair / median divergence |')
    w('|---|---|' + '---|' * (len(cols) + 1))
    for t in ('ilmn', 'hifi', 'ont'):
        for r in REFS:
            x = m.loc[(t, r)]
            last = pct(x.frac_proper_pair) if t == 'ilmn' else pct(x.gc_divergence_median)
            w('| %s | %s | ' % (TNAME[t], RNAME[r]) + ' | '.join(pct(x[c], d) for c, _, d in cols) + ' | %s |' % last)
    x = m.loc[('ilmn', 'NYGC_bwamem_GRCh38alt')]
    w('| Illumina | NYGC production bwa-mem 0.7.15 (GRCh38+alt+decoy+HLA), context only | %s | %s | %s | %s | n/a | %s | %s |' % (
        pct(x.frac_unmapped), pct(x.frac_mapq0), pct(x.frac_mapq_lt20), pct(x.frac_clipped), pct(x.edit_rate_pooled, 3), pct(x.frac_proper_pair)))
    w('')
    w('More complete references expose real paralogy that GRCh38 hides: relative to GRCh38, the oracle promotes %s Illumina reads from '
      'MAPQ<20 to >=20 but demotes %s (%.1f%% of GRCh38 MAPQ>=20 reads); CHM13 promotes %s and demotes %s, and MAPQ0 in SDs rises from %s '
      '(GRCh38) to %s (oracle). MAPQ<20 is therefore a net, not a monotone, quality signal; the oracle defines the attainable level.\n' % (
          '{:,}'.format(int(tr.loc[('ilmn', 'oracle'), 'rescued_to_ge20'])), '{:,}'.format(int(tr.loc[('ilmn', 'oracle'), 'lost_to_lt20'])),
          100 * tr.loc[('ilmn', 'oracle'), 'lost_to_lt20'] / tr.loc[('ilmn', 'oracle'), 'grch38_ge20'],
          '{:,}'.format(int(tr.loc[('ilmn', 'CHM13'), 'rescued_to_ge20'])), '{:,}'.format(int(tr.loc[('ilmn', 'CHM13'), 'lost_to_lt20'])),
          pct(st.loc[('ilmn', 'SD', 'GRCh38'), 'frac_mapq0'], 1), pct(st.loc[('ilmn', 'SD', 'oracle'), 'frac_mapq0'], 1)))
    w('### MAPQ<20 by stratum and fraction of the oracle gap recovered\n')
    w('| Tech | Stratum (reads) | GRCh38 | oracle | CHM13: gap recovered | own calls: gap recovered | 1kGP phased: gap recovered |')
    w('|---|---|---|---|---|---|---|')
    for t, strata in (('ilmn', ['all', 'all_noCEN', 'CEN', 'SD', 'TR', 'other']), ('hifi', ['all', 'all_noCEN', 'CEN', 'SD', 'nonSD']), ('ont', ['all', 'all_noCEN', 'CEN', 'SD', 'nonSD'])):
        for s in strata:
            gx, ox = C(t, s, 'GRCh38'), C(t, s, 'oracle')
            w('| %s | %s (%s) | %s | %s | %s | %s | %s |' % (TNAME[t], s, '{:,}'.format(int(gx.n_reads)), pct(gx.value), pct(ox.value),
                                                       rec(t, s, 'CHM13'), rec(t, s, 'pers2p'), rec(t, s, 'pers1KG')))
    w('')
    w('### Reference bias at clean truth het SNVs (GRCh38-REF fraction, MAPQ>=20)\n')
    w('| Tech | Stratum (sites) | GRCh38 | CHM13 | own calls | 1kGP phased | oracle |')
    w('|---|---|---|---|---|---|---|')
    for t in ('ilmn', 'hifi', 'ont'):
        for s in ('ALL', 'SD', 'TR', 'other', 'CEN'):
            w('| %s | %s (%s) | ' % (TNAME[t], 'all non-CEN' if s == 'ALL' else s, '{:,}'.format(int(AB(t, s, 'GRCh38').n_sites))) +
              ' | '.join(abtxt(t, s, r) for r in REFS) + ' |')
    w('')
    w('Illumina sites (depth>=10) with zero ALT reads at MAPQ>=20: all non-CEN %s (GRCh38), %s (CHM13), %s (own calls), %s (1kGP phased), '
      '%s (oracle); SD %s, %s, %s, %s, %s; TR %s, %s, %s, %s, %s. CHM13 carries the GRCh38-ALT allele at %s of these het sites, '
      'so part of its gain is allele choice rather than completeness.\n' % tuple(
          [pct(drop('ilmn', 'ALL', r)) for r in REFS] + [pct(drop('ilmn', 'SD', r)) for r in REFS] + [pct(drop('ilmn', 'TR', r)) for r in REFS] +
          [pct(ab[(ab.tech == 'ilmn') & (ab.ref == 'CHM13') & (ab.stratum == 'ALL') & (ab.siteset == 'clean')].frac_sites_refgenome_carries_alt.iloc[0], 1)]))
    w('### Quality of the personalization inputs (chr20 confident region minus centromere model, %s bp)\n' % '{:,}'.format(int(cs.iloc[0].confident_bp)))
    w('| Call set | het SNV recall | hom-ALT SNV recall | SNV precision | phase agreement with trio assembly | adjacent switch rate |')
    w('|---|---|---|---|---|---|')
    for k, nm in (('twopass_bcftools_ilmn', 'two-pass bcftools (same Illumina reads)'), ('1kGP_phased_panel', '1kGP phased panel')):
        x = cs.loc[k]
        pa = pct(x.phase_agree_frac, 2) if 'phase_agree_frac' in x and pd.notna(x.phase_agree_frac) else 'unphased'
        sw = pct(x.switch_rate_adjacent, 3) if 'switch_rate_adjacent' in x and pd.notna(x.switch_rate_adjacent) else 'n/a'
        w('| %s | %s | %s | %s | %s | %s |' % (nm, pct(x.recall_genotype_het, 1), pct(x.recall_genotype_hom, 1), pct(x.precision_allele, 2), pa, sw))
    rs_path = os.path.join(R, 'L2L3_callset_recall_by_stratum.tsv')
    if os.path.exists(rs_path):
        rs = pd.read_csv(rs_path, sep='\t')
        rs = rs[rs.siteset == 'clean']
        w('')
        w('| Stratum (clean truth het SNVs) | present in two-pass calls | present in 1kGP phased panel |')
        w('|---|---|---|')
        for s_ in ['other', 'TR', 'SD', 'CEN']:
            x2 = rs[(rs.stratum == s_) & (rs.callset == 'twopass_bcftools_ilmn')].iloc[0]
            x1 = rs[(rs.stratum == s_) & (rs.callset == '1kGP_phased_panel')].iloc[0]
            w('| %s (%s) | %s | %s |' % (s_, '{:,}'.format(int(x2.truth_het_snv)), pct(x2.present_frac, 1), pct(x1.present_frac, 1)))
    w('')
    w('### Long-read residual errors versus the own assembly\n')
    w('| Context | HiFi errors/kb | HiFi del events/kb | ONT errors/kb | ONT del events/kb | share of assembly bases |')
    w('|---|---|---|---|---|---|')
    for c in ['HP1', 'HP2-3', 'STR', 'HP4-5', 'HP6-7', 'HP8-10', 'HP11-15', 'HP16+', 'ALL']:
        w('| %s | %.2f | %.2f | %.1f | %.1f | %s |' % (c, hc.loc[c, 'err_per_kb'], hc.loc[c, 'del_ev_per_kb'], oc.loc[c, 'err_per_kb'],
                                                    oc.loc[c, 'del_ev_per_kb'], pct(hc.loc[c, 'frac_of_bases'], 1)))
    w('')
    w('*HP1 = assembly bases outside homopolymer runs and outside 2-6 bp-unit STRs; HPk = bases inside a homopolymer run of length k; deletion events are attributed to the first deleted assembly base.*\n')
    w('| Reported QV bin | HiFi bases | HiFi empirical QV | ONT bases | ONT empirical QV |\n|---|---|---|---|---|')
    for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 100)):
        hb = qv['hifi'][(qv['hifi'].reported_qv >= lo) & (qv['hifi'].reported_qv < hi)].bases.sum()
        ob = qv['ont'][(qv['ont'].reported_qv >= lo) & (qv['ont'].reported_qv < hi)].bases.sum()
        w('| %d-%d | %.1f%% | %.1f | %.1f%% | %.1f |' % (lo, hi - 1, 100 * hb / qv['hifi'].bases.sum(), qbin('hifi', lo, hi),
                                                     100 * ob / qv['ont'].bases.sum(), qbin('ont', lo, hi)))
    w('')
    w('![Figure 2](%s)\n' % L.get('L2L3_fig2_longread_residuals.png', 'L2L3_fig2_longread_residuals.png'))
    w('*Figure 2. Left: long-read errors per kb of the sample\'s own assembly by assembly context (homopolymer run length, 2-6 bp-unit '
      'STR), MAPQ>=20 alignments of chr20 reads to their best haplotype. Right: reported base QV versus empirical QV (mismatch and '
      'insertion errors per aligned base; bins with >=10k bases); dashed line = perfect calibration.*\n')
    txt = '\n'.join(out) + '\n'
    if a.tail:
        tail = open(a.tail).read()
        for f, target in L.items():
            tail = tail.replace('`%s`' % f, '[%s](%s)' % (f, target))
        txt += '\n' + tail
    if L:
        txt += '\n## Files\n\n' + '\n'.join('- [%s](%s) - %s' % (f, L[f], D.get(f, '')) for f in L) + '\n'
    open(a.out, 'w').write(txt)


if __name__ == '__main__':
    main()
