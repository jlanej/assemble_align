# copyproof: implementation guide

copyproof proofreads released genome assemblies for copies that are missing, doubled, joined to the wrong neighbour or placed on the wrong haplotype. It is the detection half of proposal F4 in the report ([`../report/assemble_align_review.md`](../report/assemble_align_review.md), section F4). The repair half, changing how an assembler walks its graph, is a separate and harder project that can start from copyproof's issue lists.

**Status.** Design only; nothing is implemented. This guide describes what to build, why, with which public data, and how to tell whether it works. It is meant to be steered by results: the numbers in it are starting points, and §9 lists the decisions to take once data are in.

## How to use this guide

- §1–§3 give the idea and the evidence, §8 says how success is judged, and §10 suggests an order of work.
- Treat it as a map, not a contract. When the data disagree with it, follow the data and record why in [NOTEBOOK.md](NOTEBOOK.md).
- The cheapest experiment that could stop the project, false flags on a truth genome (§8, rung 3), is placed early on purpose.
- Conventions: BED intervals are 0-based and half-open. A *haplotype* is one released haplotype FASTA together with the unplaced sequence released with it. An *issue* is a flagged locus; it becomes an *error* only after adjudication.

## 1. Why build it

### 1.1 A failure class that coverage QC cannot see

Released human assemblies are checked mainly through read coverage. HMM-Flagger models the coverage of reads aligned to the diploid assembly. NucFlag looks for secondary bases and coverage breaks in HiFi alignments. Merqury summarises k-mer QV, completeness and phasing genome-wide. All three ask whether something at a position looks wrong. A copy that is real but sits in the wrong place has normal coverage and passes.

Two 1000 Genomes cases re-tested in this repository show the class ([`../litmus/cases_HG01943_HG01103/`](../litmus/cases_HG01943_HG01103/); report chapter 15):

- **HG01943, HPRC release 2.** A de novo maternal 1.09-Mb tandem duplication (chr2:46,970,725–48,057,041) is supported by read copy number 3.05 (CRAM-index depth) to 3.13 (k-mers), 31 HiFi junction reads and 35 short-read junction pairs or splits, and by nothing in either parent. The released haplotypes contain none of the 28 junction-spanning 31-mers. Two maternal fragments were scaffolded into the paternal chromosome at four 100-N joins that no read spans. One paternal piece carries 5,638 maternal and 0 paternal inherited hapmers, and the joins created two false duplications. The paternal pieces are ordered exactly as their homologous starts on the maternal sequence, the signature of reference-guided ordering with `--dual-scaf`. The default (conservative) HMM-Flagger tracks flag 0 of 5 scored targets; the unfiltered HiFi and ONT predictions and NucFlag reach 2 of 5.
- **HG01103, 1KGP PacBio trio hifiasm.** A sub-clonal chr2/chr14 chromothripsis with 16 junctions is present in about 64% of cells. The FASTA holds all 16 junctions, but at 11 of 16 forks the long contigs follow the minority clone, and at least one contig switches clones. A review that reads the published assembly-to-GRCh38 BAM sees 5 of 16 junctions, because the BAM lacks 75 of 217 paternal contigs. The lost segments carry 0.42–0.46 of haploid HiFi depth, close to what a coverage model calls a false duplication.

Both loci were found because someone looked. How often this class occurs across a release is unknown. A tool that measures it, and a curated issue list for HPRC release 2 and HGSVC3, are both worth having.

### 1.2 What exists, and what would be new

No single signal proposed here is new:

| Tool | Evidence | Needs | Gap for this class |
|---|---|---|---|
| Merqury (Rhie 2020) | k-mer QV, completeness, hapmer phasing, genome-wide | short-read k-mers; parents for phasing | genome-wide summaries: 28 missing k-mers vanish among billions |
| Merfin K* (Formenti 2022) | expected vs observed copy number per k-mer (positive = collapsed, negative = expanded) | read k-mers and the assembly | the closest prior art for S1 below; per-k-mer values, with no segment calls and no join, piece or junction logic |
| HMM-Flagger (Asri 2026, preprint) | coverage HMM over reads aligned to the diploid assembly | HiFi or ONT alignments | published for every HPRC R2 haplotype; a misplaced real copy has normal coverage |
| NucFlag (Oshima 2025) | secondary bases and coverage; misjoin calls | HiFi alignments | dense genome-wide (HG01943: 1,054 MISJOIN calls on the paternal haplotype, about 361 per Gb, and a call at every one of the 42 gaps) |
| GAVISUNK (Dishuck 2022) | inter-SUNK distances in ONT reads | ultra-long ONT | needs ultra-long reads; weak where copies are near-identical |
| Inspector (Chen 2021), CRAQ (Li 2023) | clipping and structural signals from reads realigned to the assembly | reads | local; no copy-number, trio or release-level reasoning |
| Comparison with a truth genome (Q100 benchmark, Hansen 2026; GQC runs) | assembly against the truth | a truth genome | exact labels, but only for HG002 |

copyproof would add:

1. Segment-level calls of copy-number discordance, with an explicit "cannot assess" mask.
2. Localised recovery of read k-mers that are missing from the assembly, classified as a missing junction or a missing allele.
3. A piece and join audit: does every scaffold piece carry the right parent's hapmers, and does any read span each join?
4. One evidence card per issue, lifted to GRCh38 and CHM13.
5. Scanning of whole releases, with a curated issue list and a prevalence estimate.

Put plainly, copyproof is mostly a thin integrator. The join, piece and junction audit is the new part, and the release-scale prevalence is the main result.

### 1.3 Who would use the output

- Release curators (HPRC, HGSVC): a triaged list of issues to fix or annotate.
- Users of pangenomes and variant catalogues: a per-haplotype mask of unreliable intervals, lifted to GRCh38 and CHM13.
- Assembler developers: recurrent patterns, such as dual-scaffold joins, with worked examples.

## 2. Issue classes

Each issue is a hypothesis with evidence, and one locus can carry several classes. Treat the list as a starting vocabulary: merge or split classes when the data say so.

| Class | Meaning | Main evidence | Example | Rule out first |
|---|---|---|---|---|
| COLLAPSE | fewer copies in the assembly than in the genome | read CN > assembly CN over a segment | HG01943: 17 bins at assembly CN 2 against read CN ≈ 3 | unplaced copies not counted; GC bias; source discordance |
| EXPANSION | a false duplication | assembly CN > read CN; neighbouring pieces overlap | HG01943: overlaps A/C and C/E; 3 bins at CN 4 | a deletion in the reads' DNA; mosaic loss |
| MISSING_JUNCTION | a breakpoint present in the reads is absent; its flanks are far apart | solid read k-mers absent from the assembly; flanks map to distant loci | HG01943: 0 of 28 junction 31-mers in the release | chimeric reads; contamination |
| MISSING_SEQUENCE | an allele or insertion present in the reads is absent; its flanks are adjacent | absent k-mers with local flanks | none yet | a variant on the other haplotype; base errors |
| MISPLACED_PIECE | a scaffold piece carries the other parent's sequence | inherited hapmers per piece | HG01943 piece B: 5,638 maternal vs 0 paternal | low hapmer density (runs of homozygosity); parental sequencing error |
| UNSUPPORTED_JOIN | a gap or adjacency that no read supports | no flank-anchored spanning read; clip clusters; neighbour overlap | HG01943: four 100-N joins, none spanned | joins inside repeats longer than the reads (these become CANNOT_ASSESS) |
| SUBCLONAL | real sequence present in a fraction of cells | depth ratio stable at one fraction; absent from both parents | HG01103: 0.42–0.46 of haploid depth over 46.2 Mb | a false duplication; a copy loss |
| CLONE_SWITCH | a contig path changes clone at a fork | branch read fractions disagree along the path | HG01103: contig h1tg000031l | — |
| SOURCE_DISCORDANT | the reads' DNA and the assembly's DNA differ | discordance in NYGC Illumina but not in the assembly's own long reads | none yet | a true error (decide with the assembly's own reads) |
| CANNOT_ASSESS | the evidence cannot decide here | the mask (§5) | satellites, rDNA, high-copy k-mers | — |

## 3. Evidence channels

For each channel: the idea, what the pilot did, what is known, the open choices and the pitfalls. Pilot scripts are in [`../litmus/cases_HG01943_HG01103/`](../litmus/cases_HG01943_HG01103/).

### S1. Positional k-mer copy number

**Idea.** In windows along each released haplotype, compare the copy number of the window's k-mers in the reads with their copy number in the sample's whole release. The whole release means both haplotypes plus every unplaced, unassigned, disconnected, rDNA, EBV or mitochondrial sequence. Read CN is the k-mer's count divided by the haploid k-mer depth; assembly CN is the k-mer's count across the release. Flag runs of windows where the two disagree.

**Pilot.** `kmer_cn_track.py` reads meryl-lookup `-wig-count` tracks for the child, both parents and the assembly. Per 10-kb bin it reports two values. The first is the median read count of k-mers present in both reads and assembly, divided by haploid depth; the second is the median assembly count. On HG01943 it gives 3.13 read copies against assembly CN 3 in 87 bins, 2 in 17 bins (the missing copy) and 4 in 3 bins (the false duplications).

**Known.** Measured on a 2.25-Mb slice only. The genome-wide false-flag rate is unmeasured (§8, rung 3).

**Open choices.**
- *Read source.* NYGC Illumina reads are independent of the assembly, with a different DNA preparation and platform. They exist for every 1kGP sample, and for both parents in 602 trios. HiFi reads from the assembly's own DNA avoid source discordance, but they built the assembly. A reasonable start is Illumina for detection and long reads for adjudication.
- *k.* A larger k keeps more k-mers unique inside segmental duplications; a smaller k gives more solid k-mers at 30x. Measure the fraction of informative k-mers per window, by stratum, for k = 21, 25 and 31.
- *Statistic.* Window medians are robust. A per-window likelihood under integer copy number gives calibrated evidence and copes with sparse windows; it uses a negative-binomial count model whose overdispersion is fitted genome-wide.
- *Segmentation.* Run-length rules (as in the pilot), circular binary segmentation, or an HMM over integer copy-number differences. The HMM fits the count model naturally.
- *Haplotype attribution.* The comparison is diploid. Assigning a discordance to one haplotype needs more evidence: hapmers (S3), read phasing, or which haplotype holds the extra or missing copy.

**Pitfalls.**
- Estimate haploid depth from each read set's k-mer spectrum (the heterozygous peak), and correct for GC.
- In males, X and Y have a haploid expectation.
- High-copy k-mers saturate; mask them (§5).
- Releases differ in how they ship EBV, mitochondrial and unplaced sequence. Count all of it, or the comparison invents collapses.

**Prior art as a baseline.** Merfin's K* is the per-k-mer version of S1. Running Merfin on the same inputs is a cheap baseline. If K* segments as well as a custom model, use it.

### S2. Localised missing k-mers

**Idea.** Collect solid read k-mers that are absent from the whole release. Most are expected: sequencing errors near the solidity threshold, residual base errors, contamination and EBV. Cluster the rest by co-occurrence in reads, assemble each cluster into a short contig, and map its flanks. Flanks far apart (different contigs, distant positions or opposite orientations) indicate a missing junction. Adjacent flanks indicate a missing allele or a collapsed variant copy.

**Pilot.** `junction_from_reads.py` derives junction k-mers from reads. `asm_kmer_search.py` searches whole assemblies (Aho–Corasick), and `asm_junction_status.py` reports each junction's contig context. In HG01943 none of the 28 junction 31-mers is in the release, and a local hifiasm 0.25.0 run recovers them in a 44-kb contig. In HG01103 the trio assembly FASTA holds all 16 junctions.

**Known.** Merqury's completeness score cannot see 28 missing k-mers among billions; the localised version can. The number of clusters per genome is unknown.

**Open choices.**
- The solidity threshold, for example a count of at least half the haploid depth.
- A trio filter: a k-mer present in a parent is inherited; one absent from both parents is de novo or an artefact.
- The minimum cluster size: a breakpoint yields up to k − 1 spanning k-mers, and 22–28 per junction in the two cases.
- How to collect each cluster's reads: one streaming pass over the CRAM or FASTQ against a k-mer set, including unmapped reads.
- The local assembler.

**Pitfalls.**
- Derive junction k-mers from reads, never by concatenating reference flanks. The HG01943 junction carries a 2-bp CT insertion that reference concatenation would miss.
- N-gaps create missing k-mers by design.
- Source discordance creates genuinely missing k-mers that are not assembly errors.

### S3. Piece-level parental hapmer audit

**Idea.** Split each haplotype into pieces at N-runs, and at contig boundaries where an AGP or gap BED gives them. For each piece, count inherited maternal-only and paternal-only hapmers, that is, k-mers present in the child and in only one parent. Flag pieces dominated by the other parent's hapmers.

**Pilot.** `trio_phase.py` and `read_parent.py`, with the rule ≥ 100 hapmers at ≥ 3:1. Piece B on the HG01943 paternal haplotype carries 5,638 maternal against 0 paternal hapmers.

**Open choices.** The hapmer k (yak uses 31). Thresholds that scale with piece length and local hapmer density. Sliding windows inside contigs as well as whole pieces, which overlaps Merqury's switch-error statistics.

**Coverage.** The audit needs both parents. From [`../data/sample_overlap.tsv`](../data/sample_overlap.tsv): 121 of the 223 HPRC R2 samples in the 1kGP set, and 30 of the 63 HGSVC3 samples, have both parents in NYGC 30x. For the rest, consider long-read haplotags or methylation-guided phasing, or report phase as not assessable.

**Pitfalls.**
- Hapmer density is near zero in runs of homozygosity; report those as "not assessable", not "consistent".
- Build parental databases from whole genomes. In the L4 pilot, chr20-sliced parental k-mers gave 0.9–3.9% hamming error on the truth haplotypes themselves.
- A male's X carries no paternal hapmers.

### S4. Join audit

**Idea.** For every N-gap and contig adjacency, ask whether any read links the two flanks. A link is a read anchored by at least L kb on each side (pilot: 5 kb, MAPQ ≥ 20), with consistent orientation and distance. Self-align each piece against its neighbours: an overlap of more than a few kb means the join duplicated sequence. Record the order of pieces relative to the other haplotype, to catch the dual-scaffold ordering signature seen in HG01943.

**Pilot.** `asm_read_qc.py`. In HG01943 none of the four joins is spanned, and the overlaps A/C and C/E are the two false duplications.

**Inputs that already exist.** For every HPRC R2 sample, HPRC publishes the HiFi and ONT alignments to the diploid assembly that Flagger used (minimap2 2.28, "corrected" BAMs), and a gap BED per haplotype. Both are indexed in `human-pangenomics/hprc_intermediate_assembly/data_tables`. Reusing them avoids realignment. HGSVC3 assemblies would need their own alignments.

**Pitfalls.**
- Check that an alignment file covers every contig of the FASTA before trusting an absence of evidence; this is the HG01103 lesson.
- Reads cannot span joins inside repeats longer than the reads. Call those CANNOT_ASSESS, or bring in ultra-long ONT reads.
- In homozygous regions, reads aligned to a diploid assembly split between the haplotypes at MAPQ 0, so count across both haplotypes.

### S5. Clip clusters

**Idea.** Cluster the soft-clip positions of long reads aligned to the assembly. Clusters at both ends of a segment, together with a missing junction, are the signature of a collapsed tandem copy.

**Pilot.** In HG01943, 31 and 30 junction reads end at the two ends of the maternal copy. HPRC also publishes a per-sample high-clip bigWig from Flagger's processing.

**Pitfalls.** Clipping is common at heterozygous SVs, in centromeres and low-complexity sequence, and at ONT adapters. Treat clip clusters as supporting evidence only.

### S6. Cross-assembly and external evidence (optional)

- Published Flagger and NucFlag calls, used as features rather than as truth.
- A second assembly of the same sample:
  - HGSVC3 publishes Verkko and hifiasm assemblies built from the same reads, which isolates assembler disagreement.
  - Four samples (NA19240, HG00733, HG02818, NA19036) also have HPRC R2 and HGSVC3 assemblies built from different data.
- Long reads from other DNA preparations (1kGP-LRSC ONT, 1KGP PacBio, 1kG-ONT Vienna), to settle source discordance and estimate clone fractions.
- Recurrence across samples. The same misplacement pattern in many samples suggests assembler behaviour; the same junction in many samples suggests a common SV rather than an error.

## 4. From evidence to calls

- Start with transparent per-class rules, as the pilot did, with every threshold in one configuration file.
- Once labelled data exist (§8, rungs 2 and 4), fit a calibrated score per class. A small, interpretable model on the channel features and the stratum is enough, such as logistic regression or gradient-boosted trees. Calibrate it with isotonic regression or a conformal method, so that a score of 0.9 means about nine in ten such calls are real.
- The report's proposal U9 gives orientation targets: an area under the precision–recall curve at least 1.25 times that of the best single detector, and a calibration slope of 0.9–1.1.
- Keep the rules as a baseline. If the model does not beat them, ship the rules.

## 5. Where it cannot judge

Report CANNOT_ASSESS rather than guessing. Candidate parts of the mask:

- satellite arrays and centromeres (the HPRC centromere and satellite annotation per assembly; the HGSVC3 centromere release);
- rDNA;
- k-mers whose multiplicity in the release exceeds a cap;
- windows with too few informative k-mers;
- extreme GC, where Illumina coverage drops;
- a flank around each N-gap (for S1).

Report every metric by stratum:

- unique sequence;
- segmental duplications by identity (the HPRC segdup annotation per assembly);
- tandem repeats;
- pericentromeric sequence;
- chrX, chrY and the PARs.

The masked fraction of each genome is itself a reported result.

## 6. Inputs and outputs

**Inputs per sample.**

| Input | Typical source | Needed for |
|---|---|---|
| All released FASTAs of the sample, with `.fai` | HPRC R2 index; HGSVC3 assembly folders | everything |
| The sample's short-read k-mer database | NYGC 30x CRAM, decoded with GRCh38_full_analysis_set_plus_decoy_hla | S1, S2 |
| The parents' short-read k-mer databases | the parents' NYGC 30x CRAMs | S3; the S2 trio filter |
| Long reads aligned to the assembly | HPRC Flagger BAMs, or your own HiFi or ONT alignments | S4, S5, adjudication |
| A gap BED per haplotype | the HPRC gaps index, or N-runs | S3, S4 |
| Annotation (segdups, censat) | HPRC annotation indexes; the HGSVC3 release | mask, strata |
| Liftover | HPRC chains to GRCh38 and CHM13, or the impg all-vs-1 PAF | outputs; comparison across samples |
| Published QC calls | the HPRC Flagger and NucFlag indexes | features; comparison |

**Outputs (suggested).**

- **`issues.bed`**, one per haplotype: contig, start, end, issue id, score and class, as BED9 coloured by class for browsing. The R2 scorer (`score_flagged_cases.py qc-hg01943 --bed calls.bed`) reads contig, 0-based start, end and label, so this file can be scored as it is.
- **`issues.tsv`**: one row per issue, with every evidence field:
  - read and assembly CN, and the number of windows;
  - missing k-mers and flank positions;
  - hapmers by parent, spanning reads and clip counts;
  - overlap with Flagger and NucFlag calls, masked fraction and stratum;
  - GRCh38 and CHM13 coordinates;
  - a curation status (open, confirmed, rejected or source discordance).
- **An evidence card per issue** (JSON plus a small figure): the CN track, hapmer bars per piece, join support, the clip pileup and the missing-k-mer flanks. Cards serve curators and are also the fastest view for debugging.
- **Tracks**: read CN and assembly CN (bedGraph or bigWig), and the mask (BED).
- **`run.json`**: inputs with checksums, tool versions, parameters and timings, enough to rerun.

## 7. A suggested shape for the code

These are suggestions; profiling and data can overrule them.

- **Language.** Python orchestration (pysam, numpy, pandas) around established tools: meryl or KMC for counting, yak for hapmers, minimap2 and samtools. Move a stage to a compiled language only if profiling shows it dominates and no existing tool does the job.
- **Stages.** One subcommand per stage, each caching its outputs under a key built from input checksums:

  | Subcommand | Does |
  |---|---|
  | `prepare` | k-mer databases |
  | `cn` | S1 |
  | `missing` | S2 |
  | `pieces` | S3 |
  | `joins` | S4 |
  | `clips` | S5 |
  | `call` | rules or calibrated score |
  | `cards` | evidence cards |
  | `lift` | GRCh38 and CHM13 coordinates |
  | `scan` | fan-out over a release |
  | `score` | regression scorers and benchmark metrics |
- **Two modes.** Region mode works on a list of windows and suits development and the cases; the case drivers already run region-sliced in about ten minutes on 2 threads. Whole-genome mode is for real runs.
- **Release scans.** A workflow manager (Snakemake or Nextflow), with one job per sample.
- **Layout.** Once building starts, one layout that works inside this directory:
  - `src/copyproof/`: the package;
  - `tests/`: unit and synthetic tests;
  - `scripts/`: data fetching, error planting and scans;
  - `experiments/<date>_<name>/`: a README with the question, commands and small result tables, with large data kept outside git;
  - `docs/`: output formats and the curation protocol.
- **Data hygiene.** Reads, alignments, k-mer databases and assemblies stay out of git; the repository `.gitignore` already excludes the usual extensions. Scratch space comes from `mktemp` and is removed on exit, as in the case drivers.

## 8. Validation ladder

Each rung answers one question. Move up when the answer justifies the next rung's cost, and go back down when it does not.

**Rung 0: synthetic.** Small diploid genomes of a few Mb, with planted events and reads simulated at known depth. These are fast checks of every channel and class that can run on every change.

**Rung 1: the two cases.** The regression tests in [`../litmus/cases_HG01943_HG01103/regression_tests.md`](../litmus/cases_HG01943_HG01103/regression_tests.md) set three targets for a QC tool:

- R2: 5 of 5 targets with at most 2 off-target calls in the locus, reported alongside genome-wide calls per Gb;
- R4: at least 15 of 16 HG01103 junctions recovered from the FASTA;
- R5: at least 90% of the 46.2 Mb of lost-segment bases labelled sub-clonal, and under 10% labelled Dup or Err.

R1 and R3 test assemblers rather than QC, but their scorers are reusable. The pilot's 5/5 on R2 is a construct check, because its signals were designed after looking at this locus. Passing here is necessary but says little about other loci.

**Rung 2: planted errors in a truth genome.** Edit T2T-HG002 Q100 v1.1 to create each class at known positions. Then run copyproof with real HG002 reads: 30x NovaSeq PCR-free Illumina, Illumina reads of the parents (HG003, HG004), and optionally HiFi and ONT. The reads come from the true genome, so every planted edit should be found and nothing else should change. A starting design:

- **Classes.**
  - COLLAPSE: delete one copy of a segmental-duplication pair or tandem repeat.
  - EXPANSION: duplicate a segment in place or as an extra contig.
  - MISPLACED_PIECE: replace a segment with the homologous segment of the other haplotype, or insert it with N-gap joins as in HG01943.
  - UNSUPPORTED_JOIN: cut, reorder or invert, and rejoin with 100 N.
  - MISSING_SEQUENCE: remove a heterozygous insertion.
- **Sizes.** 10 kb, 50 kb, 200 kb and 1 Mb.
- **Strata.** Unique sequence; segmental duplications below 98% identity; segmental duplications at 98% identity or above.
- **Replicates.** About 10 per cell, spaced far apart so that one edited genome per class is enough.

Write a truth BED in edited coordinates. Measure sensitivity per cell, and the size and identity at which each channel stops working.

**Rung 3: false flags on the truth (go/no-go).** Run the unedited Q100 v1.1 with the same reads. Every issue is then one of three things: a tool false positive, a residual truth error, or a real difference between the reads' DNA and the assembly's DNA. The Q100 project publishes its known residual issues as `hg002v1.1_issues.v4.bed`. It holds 3,172 intervals covering 38.8 Mb, labelled DeepTrio, HapSwitch, DeepPolish, NucFlag, ErrorKmer, rDNA, Sniffles, GAP and Submitted; a breakdown of these issues underlies Table S6 of Hansen 2026. An issue landing on one of these intervals is expected. One that does not is either new or a false flag.

- Report issues per Gb by class and stratum, and inspect all of them at first.
- The report's litmus is that at most 2% of anchored duplication junctions are flagged. Make the denominator precise before running, for example every boundary in v1.1 between a segmental duplication and at least 5 kb of unique flank. Add a curator-facing number: issues per haplotype genome outside the mask.
- If the rate is high, the shape of the false positives says what to fix (§9).
- As a cheap first pass, run S1 and S4 on a few hundred megabases of slices rich in segmental duplications before the whole genome.
- v1.2 is public, with a chain to v1.1. Loci corrected between the two versions are natural positives: they should be flagged in v1.1 and not in v1.2.

**Rung 4: real assemblies with truth labels.** The Q100 benchmarking area (`T2T/HG002/assemblies/benchmarking/analyses/assembly_comparison/`) holds 24 HG002 assemblies compared with v1.1; 5 of the 24 are Q100 versions themselves (v0.1 to v1.1). Each folder holds the assembly FASTA, its alignments to v1.1 and the GQC outputs:

- error-type BED and VCF;
- SV BED;
- phased-scaffold and alignment-cluster BEDs;
- excluded regions.

Run copyproof on several: hifiasm, Verkko, HPRC year 1 and year 2, and the 2025 ONT-only assembly. Compare its issues with GQC's structural discrepancies, after projecting coordinates and removing excluded regions, and report precision and recall per class. This is the most realistic labelled test, and it also supplies training data for §4.

**Rung 5: agreement without truth.** Two kinds of pairs are available:

- the four samples with both an HPRC R2 and an HGSVC3 assembly;
- every HGSVC3 sample, comparing Verkko with hifiasm built from the same reads.

An issue in one assembly should usually coincide with disagreement between the two assemblies, and sequence on which they agree should rarely carry issues. Issues present in both assemblies of one sample point to shared assembler behaviour or source discordance.

**Rung 6: release scan.** The scan covers HPRC R2 (462 haplotypes; 121 1kGP samples with both parents in NYGC 30x) and HGSVC3 (65 genomes; 30 with both parents).

1. Triage by score.
2. Review a stratified sample per class blind: reviewers see evidence cards without scores.
3. Estimate prevalence per class, with confidence intervals.
4. Share confirmed issues with the consortia, as issues rather than errors, until they agree.

## 9. Decision points

| # | Decision | What to look at | Options | Current lean |
|---|---|---|---|---|
| D1 | Read source for copy number | rung 3 false flags by stratum; rung 2 sensitivity | Illumina; HiFi; Illumina plus long-read adjudication | Illumina to detect, long reads to adjudicate |
| D2 | k-mer size | informative k-mers per window by stratum; rung 2 sensitivity in SDs ≥ 98% | 21, 25, 31, or two sizes combined | 31 for junctions and hapmers; test 21–31 for CN |
| D3 | Coordinate frame | ease of curation; comparison with GQC and across samples | assembly coordinates; reference-projected | call in assembly coordinates, lift the outputs |
| D4 | Segmentation | rung 2 boundary accuracy; rung 3 false flags | run-length; CBS; HMM | an HMM with negative-binomial emissions, once the count model checks out |
| D5 | Negative-control denominator | the report's 2% litmus | as in §8, rung 3 | fix it before running rung 3 |
| D6 | Thresholds and score | rungs 2–4 | rules; calibrated model | freeze before rung 6 |
| D7 | Non-trio assemblies | rung 2 MISPLACED_PIECE sensitivity with and without parents | long-read haplotags; methylation; population panels; not assessable | decide from the measured loss |
| D8 | Reuse or build | Merfin K* vs S1; Flagger's high-clip track vs S5 | reuse; build | reuse when it performs as well |
| D9 | Cost | wall time and memory on two whole genomes | one k-mer database per read set, reused by every channel; region-mode scans | measure first; counting will probably dominate |
| D10 | Change of course | the rung 3 rate after masking and calibration | narrow the scope to scaffold joins and SD-rich regions | narrow if curators would drown |

## 10. Milestones

Each milestone ends with a short entry in NOTEBOOK.md stating what was learned and what changes.

- **M0. The starting point.** Re-run `run_hg01943_case.sh` and `run_hg01103_case.sh` and confirm that the R1–R5 baselines reproduce. Fetch and checksum the HG002 inputs, and record tool versions.
- **M1. The copy-number channel.** S1 in region mode with the count model, plus synthetic tests. Reproduce the HG01943 bins: 87 at assembly CN 3, 17 at CN 2 and 3 at CN 4, against read CN ≈ 3.1.
- **M2. Early go/no-go.** S1, and S4 if ready, on HG002 Q100 v1.1: first slices rich in segmental duplications, then the whole genome. Settle D5 and D10.
- **M3. The remaining channels.** S2–S5; R2 at 5/5 with at most 2 off-target calls in the locus; R4 and R5.
- **M4. Planted errors and calibration.** Rung 2; a first calibrated score; rung 3 genome-wide.
- **M5. Real errors.** Rung 4 on several HG002 assemblies; rung 5 on the four dual-assembly samples and a few HGSVC3 Verkko/hifiasm pairs.
- **M6. Release pilot.** 10–20 samples from HPRC R2 and HGSVC3, with evidence cards, a curation protocol tried on real issues, and the cost per sample.
- **M7. Release scan and write-up.** The full scan, prevalence estimates, a curated issue list shared with the consortia, and a methods paper. Assembler repair (F4 repair) and clone-aware assembly (F5) can start from this list.

## 11. Lessons from the pilots

- Region slices miss reads placed elsewhere, so copy number needs whole-genome k-mer counts.
- Region-sliced parental k-mer databases create a phasing noise floor: 0.9–3.9% hamming error on truth haplotypes.
- Published alignment files can be incomplete. Compare BAM headers and per-contig read counts with the FASTA index.
- Junction k-mers come from reads: the HG01943 junction carries a 2-bp CT insertion.
- NYGC CRAMs decode with GRCh38_full_analysis_set_plus_decoy_hla. CRAM indexes alone give cheap depth estimates (`crai_cn.py`).
- Use the virtual-hosted bucket name `1000g-ont.s3.amazonaws.com`. Range reads from EBI can truncate silently, so count what you fetched.
- Cell lines change in culture: HG01103's sub-clone, EBV, and NA18612, which was withheld from HPRC R2 for chr16–chr17 joins. Not every discordance is an assembly error.
- HG002 is not in the HPRC R2 QC'd set; the R2 index lists the two hg002 v1.1 haplotypes alongside GRCh38 and CHM13.
- Large intermediates stay outside git; the case drivers make `mktemp` scratch and remove it on exit.

## 12. Ideas beyond version 1

- **Graph mode.** Read the assembler's GFA, estimate unitig copy numbers and propose copy-number-consistent walks; this is the repair half of F4. In HG01943 the non-trio unitig graph already holds a 369-kb unitig at twice haploid depth.
- **Phase checks without parents.** Methylation-guided phasing (MethPhaser, HapBridge).
- **Strand-seq** (HGSVC3) to check the orientation of joins and inversions.
- **Population view.** Project issues across the 462 haplotypes with impg to find systematic hotspots.
- **Release-diff mode.** When a consortium patches an assembly, rerun only the changed contigs.
- **Interactive review.** IGV or JBrowse sessions generated from each evidence card.

## 13. Resources

[resources.tsv](resources.tsv) lists every public resource named here, with its URL, size where known, use, and the date it was checked. The most useful starting points:

- **HPRC R2 data tables** (`github.com/human-pangenomics/hprc_intermediate_assembly`, folder `data_tables/`):
  - the assembly index: 466 rows, the 462 haplotypes plus GRCh38, CHM13 and the two hg002 v1.1 haplotypes;
  - the Flagger and NucFlag indexes, including read-to-assembly BAMs and high-clip bigWigs;
  - annotation indexes for gaps, segdups, censat and RepeatMasker, and chains to GRCh38 and CHM13.
- **HGSVC3**:
  - Verkko assemblies with unassigned and disconnected sequence shipped separately;
  - hifiasm assemblies from the same reads;
  - the accessions and QV tables under `release/Assembly_Info/v2.0/`.
- **HG002**:
  - the truth versions v1.1 and v1.2, and the known-issue BED;
  - Illumina reads of HG002 and both parents;
  - 24 assemblies compared with v1.1 by GQC.
- **NYGC 30x**: the sequence indexes and the pedigree.
- **This repository**:
  - the case drivers, helpers, scorer and baselines;
  - the dataset catalogue ([`../data/data_catalog.tsv`](../data/data_catalog.tsv)) and the sample overlap matrix ([`../data/sample_overlap.tsv`](../data/sample_overlap.tsv));
  - the remote-access recipes in `../data/`.

## 14. References

- Asri et al. 2026. Evaluating genome assemblies with HMM-Flagger. *bioRxiv*. Preprint. https://doi.org/10.64898/2026.02.27.708355
- Byrska-Bishop et al. 2022. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell*. https://doi.org/10.1016/j.cell.2022.08.004
- Chen et al. 2021. Accurate long-read de novo assembly evaluation with Inspector. *Genome Biology*. https://doi.org/10.1186/s13059-021-02527-4
- Chen et al. 2025. HapBridge: A Methylation-Guided Approach for Correcting Switch Errors and Bridging Phased Blocks in Long-Read Phasing. *bioRxiv*. Preprint. https://doi.org/10.1101/2025.11.07.687303
- Cheng et al. 2021. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*. https://doi.org/10.1038/s41592-020-01056-5
- Cheng et al. 2024. Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. *Nature Methods*. https://doi.org/10.1038/s41592-024-02269-8
- Dishuck et al. 2022. GAVISUNK: genome assembly validation via inter-SUNK distances in Oxford Nanopore reads. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btac714
- Formenti et al. 2022. Merfin: improved variant filtering, assembly evaluation and polishing via k-mer validation. *Nature Methods*. https://doi.org/10.1038/s41592-022-01445-y
- Fu et al. 2024. MethPhaser: methylation-based long-read haplotype phasing of human genomes. *Nature Communications*. https://doi.org/10.1038/s41467-024-49588-0
- Hansen et al. 2026. A complete diploid human genome benchmark for personalized genomics. *Cell*. https://doi.org/10.1016/j.cell.2026.06.016
- Kokot et al. 2017. KMC 3: counting and manipulating k-mer statistics. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btx304
- Li 2018. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/bty191
- Li et al. 2023. Identification of errors in draft genome assemblies at single-nucleotide resolution for quality assessment and improvement. *Nature Communications*. https://doi.org/10.1038/s41467-023-42336-w
- Logsdon et al. 2025. Complex genetic variation in nearly complete human genomes. *Nature*. https://doi.org/10.1038/s41586-025-09140-6
- Lucas et al. 2026. HPRC2: A human pangenome reference with near-complete coverage of common genetic variation. *bioRxiv*. Preprint. https://doi.org/10.64898/2026.07.21.739710
- Nurk et al. 2022. The complete sequence of a human genome. *Science*. https://doi.org/10.1126/science.abj6987
- Oshima 2025. logsdon-lab/NucFlag: v0.3.3.1. *Zenodo* (software). https://doi.org/10.5281/zenodo.17902823
- Rautiainen et al. 2023. Telomere-to-telomere assembly of diploid chromosomes with Verkko. *Nature Biotechnology*. https://doi.org/10.1038/s41587-023-01662-6
- Rhie et al. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biology*. https://doi.org/10.1186/s13059-020-02134-9
- Schloissnig et al. 2025. Structural variation in 1,019 diverse humans based on long-read sequencing. *Nature*. https://doi.org/10.1038/s41586-025-09290-7

Software without a paper of its own: yak (`github.com/lh3/yak`), meryl (`github.com/marbl/meryl`), Merqury (`github.com/marbl/merqury`), Merfin (`github.com/arangrhie/merfin`), HMM-Flagger (`github.com/mobinasri/flagger`), NucFlag (`github.com/logsdon-lab/NucFlag`), impg (`github.com/pangenome/impg`).
