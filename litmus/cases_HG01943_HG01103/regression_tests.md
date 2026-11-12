# Regression tests from two flagged 1000 Genomes cases

Five falsifiable tests, built from HG01943 and HG01103, for assemblers, assembly-QC tools and assembly-evaluation pipelines. Each names public inputs by URL and coordinates (no data are shipped), a metric, a measured baseline, a target and a pass/fail rule.

The scoring logic is in score_flagged_cases.py. The two drivers rebuild every input from public data in a self-cleaning mktemp scratch and run the scorer:
- `run_hg01943_case.sh OUTDIR [THREADS] [WHOLE_ASM_SCAN]`;
- `run_hg01103_case.sh OUTDIR [THREADS] [ASM_SCAN] [ASM_BAMS]`.

Background and evidence are in flagged_cases.md.

## Inputs

| Input | Source (public) | Slice |
|---|---|---|
| HG01943 trio short reads | NYGC 30x CRAMs: `1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988955/HG01943.final.cram`; `.../data/ERR3242006/HG01941.final.cram`; `.../data/ERR3242007/HG01942.final.cram` (decode with GRCh38_full_analysis_set_plus_decoy_hla.fa) | chr2:46,400,001-48,650,000 |
| HG01943 HiFi | `s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/samples/HG01943/alignments/GRCh38_no_alt/hifi_reads/HG01943.GRCh38_no_alt.sorted.bam` | chr2:46,790,001-48,240,000 (7,225 primary reads) |
| HG01943 assembly | `s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01943/assemblies/release2/HG01943_{pat,mat}_hprc_r2_v1.0.1.fa.gz` (+ `.fai`, `.gzi`, `.md5`); Winnowmap BAMs under `.../release2/alignments/winnowmap/GRCh38_no_alt/` | CM089354.1:46,450,001-49,420,000; CM089362.1:46,420,001-48,680,000; JBHIIL010000075.1 |
| HG01943 published QC | HMM-Flagger v1.2.0 and NucFlag v0.3.3 BEDs under `s3-us-west-2.amazonaws.com/human-pangenomics/submissions/` (full paths in `run_hg01943_case.sh`) | whole genome |
| HG01103 trio short reads | `.../additional_698_related/data/ERR3988840/HG01103.final.cram`; `.../data/ERR3241820/HG01101.final.cram`; `.../data/ERR3241821/HG01102.final.cram` | 32 junction ends ±1 kb; CRAM indexes for copy number |
| HG01103 HiFi | `1000g-ont.s3.amazonaws.com/1KGP_PacBio_WGS/Alignment/GRCh38/HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam` | 32 junction ends ±3 kb |
| HG01103 parents ONT | `ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG0110{1,2}.hg38.cram` | 32 junction ends ±1.1 kb |
| HG01103 assemblies | `1000g-ont.s3.amazonaws.com/1KGP_PacBio_WGS/Assembly/Trio_Assembly/HG01103-W-0/HG01103-W-0.hifiasm.dip.hap{1,2}.p_ctg.gfa.fasta.gz`; `.../Assembly/Regular_Assembly/HG01103-W-0/HG01103-W-0.hap{1,2}.p_ctg.fa.gz`; `.../ASM_to_REF_BAM/trio_asm/GRCh38/HG01103-W-0.hap{1,2}.bam` (no index) | whole files, streamed |
| Event coordinates | HG01943: segment chr2:46,970,725-48,057,041. HG01103: junctions and segments as published in `github.com/jlanej/svreview-` (`examples/hg01103/junctions.bedpe`, `segments.bed`), copied as coordinates only into `data/` | - |

Use the virtual-hosted bucket name `1000g-ont.s3.amazonaws.com`; the path-style `s3.amazonaws.com/1000g-ont/...` can be refused by proxies. Range reads from EBI can truncate silently. `lr_junctions_bedpe.py` reports `NO_READS_FETCHED` for any junction where that happened; retry those.

The junction k-mer sets are read-derived, not reference concatenations:
- **HG01943** ([hg01943_junction_31mers.tsv](hg01943_junction_31mers.tsv)): 28 31-mers spanning chr2:48,057,041 + `CT` + chr2:46,970,725, present in the child's reads and absent from both parents and from GRCh38.
- **HG01103** ([hg01103_junction_31mers.tsv](hg01103_junction_31mers.tsv)): 410 31-mers, 22-26 per junction, from the majority read string at each break; the 2 that occur in GRCh38 chr2/chr14 are removed.

## Tests

### R1 Assembly of the HG01943 duplication

**What it tests.** Whether an assembler reconstructs a constitutional tandem duplication longer than any read: both copies, the junction, both in the right parental haplotype, with no scaffold joins.

**Procedure.** Assemble the 7,225 recruited HiFi reads, or assemble genome-wide and extract the contigs aligning to chr2:46.4-48.65 Mb. Optional inputs are the parents' Illumina k-mers (yak or meryl), ONT ultra-long reads and Hi-C. Then score:

`score_flagged_cases.py asm-hg01943 --pat PAT.fa --mat MAT.fa --junction-kmers hg01943_junction_31mers.tsv --ref-window win.fa --hapmers-mat mat.inh.meryl --hapmers-pat pat.inh.meryl`

For a non-trio assembly, pass the haplotype carrying the junction as `--mat`.

**Criteria.** All four must hold to PASS:
1. At least 90% of the 28 junction 31-mers are on the maternal haplotype and none on the paternal.
2. Copy-number concordance, per haplotype, in the 131 10-kb bins of chr2:46.85-48.18 Mb is ≥ 0.95. Expected copies are 2 inside the segment on the maternal haplotype and 1 everywhere else. Scaffolds are split at N-runs, and alignments are minimap2 asm5, ≥ 10 kb.
3. There are no ≥ 10-N gaps in contigs aligned to the window.
4. No scaffold piece carries the other parent's inherited hapmers (≥ 100 and ≥ 3:1).

| Assembly | Junction 31-mers mat / pat | CN concordance mat / pat | Gaps | Trio-inconsistent pieces | Verdict |
|---|---|---|---|---|---|
| HPRC release 2 (hifiasm 0.19.9, `--ul --dual-scaf`) | 0 / 0 | 0.359 / 0.496 | 4 | 1 (B: 5,638 maternal vs 0 paternal hapmers) | FAIL, 0/4 criteria |
| hifiasm 0.25.0 trio, local reads | 28 / 0 | 0.260 / 1.000 | 0 | 0 | FAIL, 3/4 |
| hifiasm 0.25.0 non-trio, local reads | 28 / 0 | 0.206 / 1.000 | 0 | 0 | FAIL, 3/4 |
| Target | ≥ 26 / 0 | ≥ 0.95 / ≥ 0.95 | 0 | 0 | PASS |

**Headroom.** Maternal concordance can rise from 0.26 to 1.0, about 97 of 131 bins; that requires placing the second 1.09-Mb copy. The local unitig graph already holds the copy-number signal: a 369-kb unitig at twice haploid depth, entered both from the maternal flank and from the 44-kb junction unitig.

### R2 QC of the HG01943 locus

**What it tests.** Whether a QC tool flags the five defects of this locus, in assembly coordinates, without flooding the locus with other calls.

**Procedure.** Run the tool genome-wide on the public assembly with the reads it needs. Then score its BED (contig, 0-based start, end, label):

`score_flagged_cases.py qc-hg01943 --bed calls.bed`

**Targets**, in HPRC release-2 contig coordinates:
- **T1, unsupported joins.** A call with a non-gap label within 2 kb of each of the four 100-N gaps on CM089354.1 (47,396,432; 47,765,591; 48,221,142; 48,522,661). Required: all four.
- **T2, trio-inconsistent piece.** At least half of piece B (CM089354.1:47,396,532-47,765,590) covered.
- **T3, missing copy.** At least 5 kb of calls on each stretch where the second copy's reads pile up: CM089362.1:47,368,148-47,454,263, and CM089354.1:48,102,850-48,201,244 or CM089362.1:47,684,402-47,782,807.
- **T4, missing junction.** A call within 5 kb of either end of M's single copy (CM089362.1:46,997,979 or 48,085,387).
- **T5, false duplications.** At least 25% of either copy of each overlap covered: A/C (CM089354.1:47,375,929-47,396,431 or 47,765,691-47,786,203) and C/E (48,202,383-48,221,142 or 48,522,761-48,541,520).

Off-target calls are those in the locus more than 5 kb from any target; they are counted and reported. Genome-wide call density should be reported alongside.

| QC output (published by HPRC unless noted) | T1 | T2 | T3 | T4 | T5 | Targets | Off-target in locus | Verdict |
|---|---|---|---|---|---|---|---|---|
| HMM-Flagger v1.2.0 HiFi, conservative (default) | 0/4 | 0/1 | 0/2 | 0/2 | 0/2 | 0/5 | 0 | FAIL |
| HMM-Flagger v1.2.0 HiFi, full prediction | 2/4 | 0/1 | 2/2 | 0/2 | 2/2 | 2/5 | 0 | FAIL |
| HMM-Flagger v1.2.0 ONT R9.4.1, conservative | 0/4 | 0/1 | 0/2 | 0/2 | 0/2 | 0/5 | 0 | FAIL |
| HMM-Flagger v1.2.0 ONT R9.4.1, full prediction | 3/4 | 1/1 | 2/2 | 0/2 | 1/2 | 2/5 | 0 | FAIL |
| NucFlag v0.3.3 HiFi | 4/4 | 0/1 | 1/2 | 0/2 | 2/2 | 2/5 | 2 calls (2.5 kb) | FAIL |
| This work (join spanning reads, piece hapmers, k-mer CN track, clip clusters on the assembly) | 4/4 | 1/1 | 2/2 | 2/2 | 2/2 | 5/5 | 0 | PASS (construct check, not independent) |

NucFlag's genome-wide burden is 1,054 MISJOIN calls on the paternal haplotype (about 361 per Gb) and 1,220 on the maternal one, with a call at every one of the 42 gaps.

**Pass rule.** 5/5 targets and ≤ 2 off-target calls in the locus. Also report genome-wide calls per Gb. The NucFlag density above is a baseline for that report, not a threshold: no specificity target is set, because the genome-wide behaviour of the new signals was not measured.

### R3 Assembly of the HG01103 rearrangement

**What it tests.** Whether an assembler represents a cell-line rearrangement present in about 64% of cells: all 16 junctions in one haplotype, threaded through long contigs, rather than one clone chosen independently at each fork.

**Procedure.** Assemble the child's HiFi, genome-wide or from reads within 100 kb of the 32 junction ends. Then score:

`score_flagged_cases.py asm-hg01103 --hap1 H1.fa --hap2 H2.fa --junction-kmers hg01103_junction_31mers.tsv`

**Criteria.**
1. All 16 junctions (≥ 80% of each junction's 31-mers) are on one haplotype and none are on the other.
2. At least 15 of the 16 are internal: on contigs of ≥ 100 kb and ≥ 100 kb from a contig end.

A third criterion needs graph output, so it is scored manually. Every fork must carry both branches, with read-support fractions consistent with f = 0.64 ± 0.10 (short reads), and no contig may switch clones.

| Assembly | Junctions hap1 / hap2 | Internal | Context of the other junctions | Verdict |
|---|---|---|---|---|
| 1KGP-PacBio trio hifiasm 0.25.0 | 16 / 0 | 5 (J6-J10) | 4 on 42-69-kb contigs (J1, J2, J3, J12); 7 within 18-64 kb of a contig end (J4, J5, J11, J13-J16); ≥ 1 cross-clone contig (h1tg000031l) | FAIL |
| 1KGP-PacBio non-trio hifiasm 0.25.0 | 6 / 7 | not scored | J4, J5 and J7 absent | FAIL |
| Target | 16 / 0 | ≥ 15 | forks annotated with branch fractions | PASS |

### R4 Alignment-only review of HG01103

**What it tests.** Whether an assembly-evaluation workflow sees what the assembly contains. An SV or structure review reading assembly-to-reference alignments must recover the junctions present in the FASTA.

**Metric.** Junctions recovered of 16, against the FASTA k-mer truth of 16/16 on hap1.

**Baseline.** The published `ASM_to_REF_BAM` gives 5/16. It holds records for 142 of 217 paternal contigs; 75 contigs (22.9 Mb) are absent, and the 11 junction-side pieces of 17.8-63.7 kb have no record.

**Pass.** At least 15/16. Re-aligning the junction contigs with minimap2 2.31 asm5 recovers them at MAPQ 60, so this test asks the pipeline to audit sequence, not just alignments.

### R5 QC of a sub-clonal loss

**What it tests.** How a coverage-based QC labels real sequence present in a minority of cells.

**Expected.** Over the nine lost segments (46.2 Mb), the paternal contigs of the trio assembly should receive 0.42-0.46 of haploid HiFi depth. That is the measured ratio of paternal to maternal haplotagged depth on GRCh38; kept segments read 0.95-0.97.

A correct QC labels these contigs as reduced cell fraction (sub-clonal, f ≈ 0.55-0.65), not as false duplication or error. It also flags the switch point of cross-clone contigs such as h1tg000031l.

**Baseline.** Not measured: no published QC exists for this assembly, and HMM-Flagger was not run here. The projected HMM-Flagger state is Dup, because 0.42-0.46 of haploid depth is close to the half-coverage it models as false duplication.

**Pass.** At least 90% of lost-segment bases are labelled sub-clonal, and under 10% are labelled Dup or Err.

## Files

Drivers. Both were validated end to end on 2026-11-11 at 2 threads with the default options: HG01943 in 9 min 48 s, HG01103 without the optional assembly scan or BAM download in 8 min 2 s. Both reproduced the values above.
- [run_hg01943_case.sh](run_hg01943_case.sh): HG01943 end to end. It builds its own read-derived set of 24 junction 31-mers (≥ 3 bases each side), against the 28 in the shipped table; criterion 1 of R1 applies to whichever set is used.
- [run_hg01103_case.sh](run_hg01103_case.sh): HG01103 end to end.

Scorer and helpers:
- [score_flagged_cases.py](score_flagged_cases.py): the R1, R2 and R3 scorer.
- [crai_cn.py](crai_cn.py), [kmer_cn_track.py](kmer_cn_track.py), [hifi_junction.py](hifi_junction.py), [sr_junction.py](sr_junction.py), [sr_junctions_bedpe.py](sr_junctions_bedpe.py), [lr_junctions_bedpe.py](lr_junctions_bedpe.py), [junction_from_reads.py](junction_from_reads.py), [asm_kmer_search.py](asm_kmer_search.py), [trio_phase.py](trio_phase.py), [read_parent.py](read_parent.py), [asm_read_qc.py](asm_read_qc.py), [asm_junction_status.py](asm_junction_status.py).

Coordinates (driver `data/` folder):
- [hg01943_junction.bedpe](hg01943_junction.bedpe);
- [hg01103_junctions.bedpe](hg01103_junctions.bedpe) and [hg01103_segments.bed](hg01103_segments.bed), coordinates only, from the jlanej/svreview- repository.

Baseline score files:
- R1: [score_R1_hprc_r2.json](score_R1_hprc_r2.json), [score_R1_local_hifiasm_trio.json](score_R1_local_hifiasm_trio.json), [score_R1_local_hifiasm_nontrio.json](score_R1_local_hifiasm_nontrio.json);
- R2: [score_R2_flagger_hifi_conservative.json](score_R2_flagger_hifi_conservative.json), [score_R2_flagger_hifi_full.json](score_R2_flagger_hifi_full.json), [score_R2_flagger_ont_conservative.json](score_R2_flagger_ont_conservative.json), [score_R2_flagger_ont_full.json](score_R2_flagger_ont_full.json), [score_R2_nucflag_v0.3.3_hifi.json](score_R2_nucflag_v0.3.3_hifi.json), [score_R2_this_work.json](score_R2_this_work.json), with [hg01943_this_work_qc_calls.bed](hg01943_this_work_qc_calls.bed);
- R3: [score_R3_test_selected_contigs.json](score_R3_test_selected_contigs.json). This run used the extracted paternal contigs that carry the junctions; the genome-wide trio and non-trio counts in the R3 table come from the full-FASTA scan behind [hg01103_junctions.tsv](hg01103_junctions.tsv).

Tool versions used here: samtools 1.25, bcftools 1.25, minimap2 2.31-r1302, hifiasm 0.25.0-r726, meryl 1.4.2, Merqury (bioconda), yak 0.1-r69, pysam 0.24.1, pyahocorasick 2.3.1.
