# Open benchmark protocol for alignment and assembly on 1kGP-linked public data

*Status: 2026-11-11. All data are public and addressed by HTTPS URL (see data_catalog.tsv). Commands marked "tested" were executed in this track on macOS arm64 with samtools/bcftools/htslib 1.24 and minimap2 2.31, one CPU thread, streaming only the needed slices.*

The protocol scores aligners and assemblers on diverse 1000 Genomes samples against truth that does not come from the method under test. Four public resources make that possible in 2026: assembly-based benchmarks for HG002 ([Hansen 2026](https://doi.org/10.1016/j.cell.2026.06.016)) that calibrate any pipeline; pedigree-validated truth for the 1kGP sample NA12878 ([Kronenberg 2025](https://doi.org/10.1038/s41592-025-02750-y)); independent haplotype-resolved assemblies of 1kGP samples from HPRC release 2 ([Lucas 2026](https://doi.org/10.64898/2026.07.21.739710), preprint) and HGSVC3 ([Logsdon 2025](https://doi.org/10.1038/s41586-025-09140-6)); and NYGC 30x Illumina data on 602 complete trios ([Byrska-Bishop 2022](https://doi.org/10.1016/j.cell.2022.08.004)) for k-mer-based, reference-free assembly checks. Every test below states its inputs, truth, metric, current baseline (measured here where possible) and a pass/fail rule, and every metric is reported overall and in GIAB v3.6 strata ([Dwarshuis 2024](https://doi.org/10.1038/s41467-024-53260-y)).

## 1. Truth hierarchy

| Tier | Truth | Samples | Use | Caveat |
|---|---|---|---|---|
| T0 calibration | GIAB HG002 v5.0q small (smvar) and structural (stvar) benchmarks from the T2T-HG002 Q100 v1.1 assembly; GRCh37, GRCh38, CHM13v2.0 | HG002 (not 1kGP) | calibrate pipelines before applying them to 1kGP samples | smvar BED covers 2,739.4 Mb of GRCh38 and 2,767.9 Mb of CHM13 (computed here); complex SVs included, so use truvari --refine |
| T1 pedigree | Platinum Pedigree truth v1.2 (small variants, SVs, tandem repeats) + GIAB HG001 v4.2.1 | NA12878 (1kGP CEU) | only 1kGP sample with curated truth | use the held-out technical-replicate HiFi BAM for unbiased tests |
| T2 assembly-derived | dipcall/PAV calls from HPRC R2 and HGSVC3 haplotypes; regions where the two independent assemblies agree | NA19240, HG00733, HG02818, NA19036 (both consortia); 223 other HPRC R2 and 59 other HGSVC3 1kGP samples (single assembly) | truth for reads, variants and assemblies on diverse genomes | assembly errors propagate; restrict to Flagger/NucFlag-clean, concordant regions |
| T3 inheritance | NYGC 30x parental Illumina k-mers (602 complete trios) | trio children | reference-free switch/hamming error (yak), hap-mer completeness (Merqury) | k-mer methods are blind inside long exact repeats |

## 2. Region-sliced access recipes (tested)

Every test runs on slices streamed with htslib over HTTPS, so a laptop never holds a whole genome of reads.

```bash
REF=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# NYGC CRAM decoded with the remote uncompressed reference (tested: NA19240 chr20:40.0-40.1 Mb, mean depth 30.97x, 29 s)
samtools depth -a -r chr20:40000000-40100000 --reference $REF https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram
# Long-read BAMs on S3 (tested: HG01103 HiFi 41.2x in 6 s; HG01258 ONT R10 41.0x in 26 s including a 51.7 MB .bai download)
samtools depth -a -r chr20:40000000-40100000 https://1000g-ont.s3.amazonaws.com/1KGP_PacBio_WGS/Alignment/GRCh38/HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam
# bgzip assemblies: download .fai/.gzi once, then random access (tested on HPRC R2 NA19240 maternal haplotype)
A=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/assemblies/release2/NA19240_mat_hprc_r2_v1.0.1.fa.gz
curl -sO $A.fai; curl -sO $A.gzi; samtools faidx --fai-idx $(basename $A).fai --gzi-idx $(basename $A).gzi $A 'NA19240#2#CM099613.1:1000001-1000060'
# tabix VCF over HTTPS (tested: 528 GIAB v5.0q smvar records in chr20:30.0-30.1 Mb)
bcftools view -H -r chr20:30000000-30100000 https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz | wc -l
```

Practical rules: run from a `mktemp -d` scratch directory with `trap 'rm -rf "$d"' EXIT`, because htslib writes downloaded indexes to the working directory; reuse that directory across calls of one session to avoid re-downloading large ONT indexes; avoid pericentromeric probe windows (chr20:30.0-30.1 Mb gave 13.3x for the same NA19240 CRAM that gives 31.0x at chr20:40 Mb); pass `-T`/`--reference` explicitly for CRAM; Vienna ONT CRAMs need their own references (1KG_ONT_VIENNA/reference/).

## 3. Alignment tests

### A1. Mapping accuracy versus MAPQ on simulated reads with known origin

Build the two NA19240 haplotypes of a test chromosome (chr20 first, then an SD-rich set such as chr1q21, chr15q11-13 and chr16p11) by applying the 1kGP phased panel genotypes with `bcftools consensus -H 1`/`-H 2 -c chain`, so read origins can be lifted back to GRCh38. Simulate 2×150 reads with wgsim and HiFi/ONT reads with PBSIM3 ([Ono 2022](https://doi.org/10.1093/nargab/lqac092)) or Badread ([Wick 2019](https://doi.org/10.21105/joss.01316)), encode true GRCh38 coordinates in read names, align with each aligner and score with `paftools.js mapeval` ([Li 2018](https://doi.org/10.1093/bioinformatics/bty191)). Report the cumulative fraction mapped and the error rate at each MAPQ threshold, plus the fraction of reads mapped wrongly at MAPQ ≥ 30 and ≥ 60 overall and inside GIAB SegmentalDuplications and LowMappability strata. Simulation captures placement error but not the sample's structural variation, which is why A2 exists.

### A2. Real reads against the sample's own haplotypes (oracle)

Align real reads twice: to GRCh38/CHM13 (baseline) and to the sample's own diploid assembly (oracle: HPRC R2 maternal + paternal haplotypes, or HGSVC3). The oracle alignment gives a reference-bias-free read error profile and, through assembly-to-reference projection (impg/wfmash PAFs in `pangenomes/freeze/release2/impg/pafs/all-vs-1/`), a "true" reference locus for each read. Headroom = misplacement rate (reads whose GRCh38 locus disagrees with the projected locus) and lost reads (unmapped or MAPQ 0 on GRCh38 but uniquely placed on the own assembly), per stratum.

### A3. Reference bias at heterozygous sites

Take heterozygous SNVs called from the sample's assemblies (dipcall; [Li 2018](https://doi.org/10.1038/s41592-018-0054-7)) inside regions where both assemblies are confident. From each alignment (GRCh38 linear, CHM13 linear, HPRC v2.1 graph via Giraffe ([Sirén 2021](https://doi.org/10.1126/science.abg8871); [Sirén 2024](https://doi.org/10.1038/s41592-024-02407-2))), compute reference-allele fraction with `bcftools mpileup -a AD` or biastools ([Lin 2024](https://doi.org/10.1186/s13059-024-03240-8)). Metrics: mean REF fraction minus 0.5, fraction of het sites with REF fraction > 0.8 or < 0.2, stratified by distance to indels/SVs and by SD membership. Graph tests on NA19240 must use the MC v2.1 `-eval` graphs, which hold out HG002, HG005 and NA19240; the standard graphs contain NA19240's own haplotypes and would leak truth.

### A4. Variant accuracy downstream of alignment

Small variants: hap.py with the vcfeval engine ([Krusche 2019](https://doi.org/10.1038/s41587-019-0054-x); [Cleary 2015](https://doi.org/10.1101/023754), preprint) against HG002 v5.0q (with `--gender male`) and against NA12878 Platinum v1.2; vcfdist for phased small plus structural comparisons ([Dunn 2024](https://doi.org/10.1186/s13059-024-03394-5)). SVs: `truvari bench --refine` ([English 2022](https://doi.org/10.1186/s13059-022-02840-6)) with the GIAB README parameters (`--pick ac --passonly -r 2000 -C 5000`, ALT=* removed for plain bench). For 1kGP samples use T2 assembly truth restricted to agreement regions. Report F1, precision and recall overall and in LowComplexity (homopolymers ≥ 7 and tandem repeats), SegmentalDuplications, MHC/KIR and XY strata.

### A5. Runtime and memory

Measure wall time, CPU time and peak RSS with `/usr/bin/time -l` (macOS) or `-v` (Linux) at a fixed thread count, separately for index construction and mapping, on the same slice set; report reads or bases per CPU-second.

## 4. Assembly tests

Whole-genome human assembly does not fit a 64 GB laptop, so assembly tests use (i) published assemblies (HPRC R2, HGSVC3, 1kG-ONT, 1KGP PacBio) evaluated remotely, and (ii) laptop-scale re-assemblies of read subsets recruited by alignment to one chromosome or a difficult locus.

| Test | Tool / method | Metric | Truth or reference |
|---|---|---|---|
| B1 base accuracy and completeness | Merqury ([Rhie 2020](https://doi.org/10.1186/s13059-020-02134-9)) with NYGC Illumina k-mers; `yak qv` | QV, k-mer completeness | sample's own 30x Illumina |
| B2 phasing | `yak trioeval` with parental NYGC 30x k-mers ([Cheng 2021](https://doi.org/10.1038/s41592-020-01056-5)); Merqury hap-mer blocks | switch and hamming error | parents (602 trios) |
| B3 structural reliability | Flagger (HiFi read-depth; [Liao 2023](https://doi.org/10.1038/s41586-023-05896-x)), NucFlag (github.com/logsdon-lab/NucFlag), GAVISUNK ([Dishuck 2022](https://doi.org/10.1093/bioinformatics/btac714)), Inspector ([Chen 2021](https://doi.org/10.1186/s13059-021-02527-4)), VerityMap ([Bzikadze 2022](https://doi.org/10.1101/gr.276871.122)) | flagged Mb (erroneous, collapsed, duplicated, haploid) | reads only |
| B4 contiguity | fai lengths; QUAST-LG ([Mikheenko 2018](https://doi.org/10.1093/bioinformatics/bty266)) or `paftools.js asmstat` | NG50, auN, NGA50 against CHM13 | CHM13v2.0 |
| B5 T2T completeness | `t2t_contigs_remote.sh` (this track) | contigs with telomeric repeats at both ends (≥50% TTAGGG/CCCTAA hexamers in 2 kb end windows; ≥1 Mb contigs) | positive control T2T-HG002 v1.1: 46/46 |
| B6 gene completeness | compleasm ([Huang 2023](https://doi.org/10.1093/bioinformatics/btad595)); `paftools.js asmgene` | % complete single-copy genes; missing multi-copy genes | primate/mammal lineage; GENCODE on CHM13 |
| B7 region-stratified accuracy | dipcall of test assembly vs truth assembly; Merqury per-region QV | QV and discordance per Mb inside SDs (HGSVC3 SD tracks), complete centromeres (HGSVC3 1,246 validated), acrocentric short arms | HGSVC3/HPRC/T2T-HG002 assemblies |
| B8 cross-assembly concordance | dipcall/PAV of HPRC R2 and HGSVC3 haplotypes of the same sample, compared with hap.py/truvari | discordant variants per Mb; Mb of disagreement | defines T2 agreement regions |

### Measured baseline: T2T contigs in HPRC release-2 haplotypes (this track)

The script counted contigs ≥ 1 Mb with telomeric repeats at both ends, after calibration on T2T-HG002 v1.1 (46/46 chromosomes detected with 1-2 kb windows, minimum observed end fraction 0.62). A haplotype has at most 23 T2T chromosomes; acrocentric short arms rarely close.

| Assembly (HPRC R2) | contigs (all) | assembly bp | contigs ≥1 Mb | T2T contigs |
|---|---|---|---|---|
| NA19240_pat_hprc_r2_v1.0.1 | 64 | 3.039 Gb | 33 | 13 |
| NA19240_mat_hprc_r2_v1.0.1 | 65 | 3.025 Gb | 33 | 15 |
| HG01943_pat_hprc_r2_v1.0.1 | 117 | 2.921 Gb | 46 | 8 |
| HG01943_mat_hprc_r2_v1.0.1 | 127 | 3.063 Gb | 48 | 8 |
| HG00733_pat_hprc_r2_v1.0.1 | 90 | 3.036 Gb | 34 | 11 |
| HG00733_mat_hprc_r2_v1.0.1 | 98 | 3.042 Gb | 33 | 13 |
| HG02818_pat_hprc_r2_v1.0.1 | 76 | 3.019 Gb | 36 | 12 |
| HG02818_mat_hprc_r2_v1.0.1 | 112 | 3.048 Gb | 44 | 11 |
| NA19036_hap1_hprc_r2_v1.0.1 | 71 | 3.049 Gb | 33 | 15 |
| NA19036_hap2_hprc_r2_v1.0.1 | 70 | 3.028 Gb | 33 | 13 |
| HG01258_pat_hprc_r2_v1.0.1 | 67 | 2.917 Gb | 30 | 14 |
| HG01258_mat_hprc_r2_v1.0.1 | 92 | 3.040 Gb | 38 | 12 |
| HG02615_pat_hprc_r2_v1.0.1 | 72 | 3.031 Gb | 34 | 10 |
| HG02615_mat_hprc_r2_v1.0.1 | 75 | 3.032 Gb | 37 | 11 |

## 5. Litmus tests and pass/fail rules

| ID | Question | Inputs (samples) | Baseline | Oracle / target | Pass rule |
|---|---|---|---|---|---|
| L1 | Does an aligner change reduce wrong placements at high MAPQ? | A1 simulated chr20 + SD loci from NA19240 haplotypes | current aligner fraction wrong at MAPQ ≥ 60 (to be measured by pilot) | 0 wrong at MAPQ ≥ 60 | ≥ 30% fewer wrong placements at equal or higher mapped fraction, in ≥ 2 of 3 replicate simulations |
| L2 | How much reference bias remains? | NA19240 NYGC 30x; truth hets from HPRC R2 + HGSVC3 agreement | GRCh38 linear bias (pilot) | alignment to own haplotypes (A2 oracle) | headroom closed by ≥ 50% in SD and indel-proximal strata, without loss in unique sequence |
| L3 | Does a pipeline transfer from HG002 to diverse 1kGP genomes? | HG002 v5.0q; NA12878 Platinum v1.2; NA19240 T2 truth | F1 on HG002 | F1 on 1kGP samples | F1 drop from HG002 to NA19240/NA12878 ≤ 0.5 percentage points per stratum |
| L4 | Do homopolymer errors limit a read type? | same-sample panel HG01258 (Illumina 30x, TruPath, ONT R10, HiFi) | measured error events/kb in HP ≥ 8 (technology_error_profiles.md) | own-haplotype alignment | reported per technology with 95% CI over ≥ 10 windows |
| L5 | Does an assembler improvement raise completeness? | NA19240 (and the 10 recommended samples) | T2T contigs above (e.g. NA19240 13/33 pat, 15/33 mat) | 23 per haplotype | ≥ 2 additional T2T contigs per haplotype with no increase in Flagger-flagged Mb and no QV loss > 1 |
| L6 | Is a truth set trustworthy where two assemblies disagree? | NA19240, HG00733, HG02818, NA19036 | B8 discordance per Mb | 0 | report and exclude disagreement regions; any benchmark claim must hold inside agreement regions |

## 6. Reporting template

Every result row records: sample, data URL and md5/ETag, slice coordinates, reference build and URL, truth set version, stratification version (v3.6), tool name, version and exact command, thread count, wall/CPU/RSS, and the metric with its denominator. Results that rest on single-window probes are labelled provisional.

## 7. Known pitfalls

- NYGC CRAMs decode only with GRCh38_full_analysis_set_plus_decoy_hla.fa; NA19240 is in the 698-related set (ERR3989410), its parents in the 2,504 set.
- 1kG-ONT sample folders use GM identifiers for NA samples; R9.4.1 POD5 can be basecalled only with Dorado ≤ 0.9.6.
- HGSVC3: the original Verkko assembly of HG00514 is corrupted (use working/20241001_verkko_HG00514_fix); the assembly table lists NA19239's continental group as "YRI".
- HPRC MC graphs: HG00272 is held out of all v2.1 graphs (suspected chrX misassembly); benchmark graphs (-eval) additionally hold out HG002, HG005 and NA19240.
- GIAB v5.0q: hap.py needs `--gender male` for HG002; truvari bench mis-handles ALT=* (filter it) and truvari refine rewrites representations.
- GIAB stratifications exist for GRCh37, GRCh38, CHM13 and HG002 haplotypes only; sample-specific strata for 1kGP genomes must come from HGSVC3 SD and centromere files or be computed from the assemblies.
- Coverage is computed differently across sources (NYGC raw bases including duplicates, HGSVC3 per 3 Gbp, HPRC per its own pipeline, DRAGEN deduplicated); compare within a source only.

## 8. Tool versions

Latest GitHub releases queried on 2026-11-11 (api.github.com); pin these or record the commit actually used. This track ran samtools/bcftools/htslib 1.24 and minimap2 2.31-r1302 from the shared `asmalign` environment.

| Tool | GitHub repository | Latest release (date) |
|---|---|---|
| minimap2 | github.com/lh3/minimap2 | v2.31 (2026-05-19) |
| bwa-mem2 | github.com/bwa-mem2/bwa-mem2 | v2.3 (2025-06-30) |
| vg (Giraffe) | github.com/vgteam/vg | v1.78.0 (2026-10-19) |
| minigraph | github.com/lh3/minigraph | v0.21 (2024-04-21) |
| Minigraph-Cactus (cactus) | github.com/ComparativeGenomicsToolkit/cactus | v4.0.0 (2026-11-09) |
| PGGB | github.com/pangenome/pggb | v0.7.4 (2025-04-23) |
| wfmash | github.com/waveygang/wfmash | v0.24.2 (2026-01-30) |
| impg | github.com/pangenome/impg | v0.5.0 (2026-07-26) |
| pbmm2 | github.com/PacificBiosciences/pbmm2 | v26.2.0 (2026-06-30) |
| hifiasm | github.com/chhylp123/hifiasm | 0.25.0 (2025-03-18) |
| Verkko | github.com/marbl/verkko | v2.3.2 (2026-03-11) |
| samtools | github.com/samtools/samtools | 1.25 (2026-10-27) |
| hap.py | github.com/Illumina/hap.py | v0.3.15 (2021-12-04) |
| RTG Tools (vcfeval) | github.com/RealTimeGenomics/rtg-tools | 3.13 (2025-05-28) |
| truvari | github.com/ACEnglish/truvari | v5.5.0 (2026-09-19) |
| vcfdist | github.com/TimD1/vcfdist | tag:v2.6.4 |
| aardvark | github.com/PacificBiosciences/aardvark | v1.1.0 (2026-10-30) |
| dipcall | github.com/lh3/dipcall | v0.3 (2021-06-06) |
| Merqury | github.com/marbl/merqury | v1.4.1 (2026-07-06) |
| yak | github.com/lh3/yak | v0.1 (2020-07-17) |
| Flagger | github.com/mobinasri/flagger | v1.2.1 (2026-11-09) |
| NucFlag | github.com/logsdon-lab/NucFlag | v1.0.2 (2026-10-29) |
| Inspector | github.com/Maggi-Chen/Inspector | v1.3 (2024-08-25) |
| compleasm | github.com/huangnengCSU/compleasm | v0.2.10 (2026-10-07) |
| QUAST | github.com/ablab/quast | quast_5.3.0 (2024-11-10) |
| biastools | github.com/maojanlin/biastools | v0.3.3 (2026-03-12) |
| PBSIM3 | github.com/yukiteruono/pbsim3 | v3.0.5 (2025-02-18) |
| Badread | github.com/rrwick/Badread | v0.4.2 (2026-04-22) |
| DeepVariant | github.com/google/deepvariant | v1.10.0 (2026-03-05) |
| Clair3 | github.com/HKU-BAL/Clair3 | v2.0.5 (2026-11-06) |
| Sniffles2 | github.com/fritzsedlazeck/Sniffles | v2.8.1 (2026-09-10) |
| Dorado | github.com/nanoporetech/dorado | v2.1.2 (2.x line; 1.3.4 maintenance) (2026-08-26) |
