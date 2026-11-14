# Rung 3: what does copyproof flag on a truth genome?

**Question.** T2T-HG002 Q100 v1.1 is as close to a correct diploid assembly as exists. Run with
HG002's own short reads, every issue is a tool false positive, a residual truth error, or a real
difference between the reads' DNA and the assembly's. How many are there, and of what kind
(GUIDE §8 rung 3; decisions D1, D5, D10)?

**Inputs.** `hg002v1.1.fasta.gz` (maternal, paternal and chrM; 6.00 Gb); HG002 NovaSeq PCR-free
reads, both files (711.6 million 151-bp reads, 107.5 Gb, 36x); the Q100 annotation
(`hg002v1.1_issues.v4.bed`, segmental duplications, centromeres, HSats). Parents were not used
here; the same genome with its parents' reads is in `../2026-11-13_hg002_trio/`.

**Run.** `scripts/run_hg002_truth.sh DATA`, then `analyse.py` here and
`scripts/sd_boundary_litmus.py`. k = 31, 1 in 4 k-mers tracked, 1-kb windows. A second run adds
`--long-read-bam` with the public alignment of HiFi reads to this assembly, read remotely.

## Result

| | |
|---|---|
| Depth per copy | 13.39 |
| Assessable | 88.7% of bases; 679 Mb masked (94.8% of centromere and HSat annotation, 38.5% of segmental duplications) |
| Depth-deviating segments found by the HMM | 2,792 |
| ... labelled depth anomaly (`DEPTH`, `IRREGULAR`), not issues | 2,630 |
| ... labelled whole-copy change (`CN`) | 20 |
| ... ambiguous (`CN_OR_DEPTH`), reported only with support | 142, of which 30 reported |
| Read k-mers absent from the assembly, solid | 9.4 million (2.6 million low-complexity, 0.4 million high-copy, 6.3 million candidates) |
| ... reads recruited | 370,563 of 711.6 million |
| ... loci dropped because every read is on one strand | 13,102 |
| ... loci kept | 1,408 |
| Runs of assembly k-mers without any read | 13,044; 1,123 with normally covered flanks |

**The litmus (denominator fixed beforehand, NOTEBOOK entry g).** Of 20,237 anchored
duplication boundaries, 23 (0.11%) have an issue within 5 kb, 18 (0.09%) a high-confidence one.
The report's target was at most 2%.

**Issues on the truth genome, short reads only** (404 rows; diploid; `partner` and `twin` rows are
the same event seen from another copy and are not counted in the last column):

| Class | High | Low | Of the high ones, on a known Q100 issue | Loci per haplotype, high |
|---|---|---|---|---|
| COLLAPSE | 2 | 0 | 0 | 1 |
| EXPANSION | 1 (+1 partner) | 38 (+2) | 0 | 0.5 |
| UNSUPPORTED | 0 | 5 | | 0 |
| BASE_ERROR | 2 | 0 | 2 | 1 |
| MISSING_HAPLOTYPE | 7 (+14 twin) | 11 (+9) | 2 | 3.5 |
| MISSING_ALLELE | 37 (+25 twin) | 55 (+45) | 4 | 18.5 |
| MISSING_JUNCTION | 0 | 1 | | 0 |
| MISSING_SEQUENCE | 0 | 43 (+106) | | 0 |
| **All** | **49** | **153** | **8** | **25** |

So the curator-facing number is about 25 high-confidence issues per haplotype genome on a truth
genome, and about 75 low-confidence ones. (An event whose flanks exist on both haplotypes is
counted once; which haplotype gets the `locus` and which the `twin` is arbitrary without other
evidence, so per-haplotype counts are averages.) Weak one-sided `MISSING_SEQUENCE` clusters, 892
of them here, are not issues; they stay in `missing_clusters.tsv`.

**With long reads** (`--long-read-bam`, the public HiFi-to-Q100 alignment; 525 rows). The audit
adds 131 issues, 126 of them on the truth genome's own list of known issues:

| Class added or changed | High | Low | On a known Q100 issue |
|---|---|---|---|
| BASE_ERROR (long reads run through, none holds the k-mers) | 61 (was 2) | 50 | 107 of 111 |
| FALSE_JOIN (no long read runs through) | 11 | 3 | 13 of 14 |
| UNSUPPORTED_JOIN (the rDNA gaps; high when another issue lies beside the gap) | 7 | 1 | 8 of 8 |
| EXPANSION, low | | 29 (was 38) | |

Ten short copy-number losses (9 `EXPANSION`, 1 `UNSUPPORTED`) are dropped because HiFi reads
carry the k-mers that short reads lack.

## What the flags are

1. **Depth anomalies are the bulk, and are not copy-number errors.** 2,630 segments, mostly 2 kb,
   where every k-mer class is depleted by one factor (2,229) or by no regular pattern (401).
   Examples: GC-rich minisatellite stretches near chromosome ends at 0.4 to 0.6 of the expected
   depth (chr2 1.52 Mb, 20 kb with one-copy k-mers at 0.55 and two-copy k-mers at 1.0); the LPA
   kringle repeat at 0.75 of expected for multiplicities 1 to 8. Dips are 60 times as common as
   excess depth (2,589 against 41). 197 anomalies lie on known Q100 issues.
2. **Whole-copy calls are few** (20 segments, `cn_candidates.tsv`).
   - chrM (hundreds of copies per cell; excluded from the issues by name).
   - The first 60 kb of chr1 on both haplotypes, where k-mers the assembly holds 3 to 7 times are
     read at 3 to 4 more copies each (`COLLAPSE`, high confidence).
   - Three stretches of acrocentric short arms and one of chr9 where k-mers of every multiplicity
     from 2 to 8 are read at exactly one copy fewer (chr13 paternal 10.04 Mb: 2, 3, 4, 5, 6, 7, 8
     copies in the assembly read as 1.4, 2.1, 3.2, 4.0, 5.5, 6.3, 7.3).
   - Thirteen stretches of 1 to 3 kb where k-mers specific to one haplotype have no short reads at
     all while the shared k-mers sit at one copy. **Long reads say what these are**
     (`hifi_at_short_losses.tsv`): HiFi reads carry every one of the unread k-mers at 10 of the
     13, so the sequence is real and short reads fail on one allele; they carry few or none at 2
     (chr5 paternal 188.49 Mb, chr18 maternal 23.78 Mb), and both of those are known Q100 issues;
     one has too few reads. This is why a whole-copy call on fewer than 5 windows is of low
     confidence without long reads, and is dropped when long reads carry its k-mers.
3. **Immunoglobulin loci** (chr14 maternal 107.6 Mb, chr2 maternal 90.4 to 90.5 Mb, chr22) show
   losses, missing haplotype stretches and the one missing junction. HG002's DNA is from a B-cell
   line, where V(D)J recombination deletes sequence: these are source discordance, the class the
   guide said to expect.
4. **Read k-mers absent from the assembly are almost all artefacts, in four groups** that had to be
   removed in turn before anything else was visible: adapter read-through; k-mers spanning short
   tandem repeats and homopolymers; single-base errors of k-mers that are read thousands of times
   (satellites, rDNA, and the mitochondrial genome, whose error k-mers alone gave 135,000 tracked
   k-mers and 228,000 recruited reads); and recurrent miscalls in fixed sequence contexts (mostly
   T to G), which come from one read strand only. 13,102 loci were dropped for strand bias against
   1,408 kept.
5. **Runs of assembly k-mers with no short read are not usable alone, and long reads sort them
   cleanly.** 1,123 runs have normally covered flanks; 948 of those are haplotype-specific
   k-mers. The question put to the HiFi reads lying over each run is whether at least two of them
   hold its k-mers (`unsupported_runs_hifi.tsv`):

   | HiFi reads ... | Runs | On a known Q100 issue | Changed in v1.2 |
   |---|---|---|---|
   | carry the k-mers | 957 | 74 | 22 |
   | do not carry them | 119 | 115 | 27 |
   | too few reads or k-mers to say | 47 | 44 | 1 |

   So 85% of the runs are correct sequence that the short-read chemistry does not read, and the
   runs that neither read type holds are, 115 times in 119, places the truth genome's curators
   already list. By alignment alone, HiFi reads "run through" 109 of the 119: an alignment across
   a place does not show that the read holds the assembly's sequence there.

   **What the reads hold instead** (`reads_allele.py`, `not_carried_reads_allele.tsv`). For each
   of the 119, the assembly's sequence between two 25-base anchors outside the unread k-mers was
   compared with what each HiFi read holds between the same anchors. In 110 runs reads hold both
   anchors (median 26 reads), and in every one of those at most one read agrees with the
   assembly; in 101, none does. The commonest version among the reads differs from the assembly by:

   | Difference | Runs |
   |---|---|
   | length of a homopolymer | 28 |
   | another insertion or deletion of 1 base | 24 |
   | insertion or deletion of 2 to 10 bases | 35 |
   | longer insertion or deletion | 19 |
   | several bases, same length | 3 |
   | one substitution | 1 |
   | no read holds both anchors (rDNA gap edges and the like) | 9 |

   These are length differences in homopolymers and short tandem repeats, where reads also
   disagree with each other: the commonest version is held by about two thirds of the reads
   (quartiles 44% and 80%). So "the assembly's version is in no read" is firm, and "the reads' commonest version is
   the right one" is not.

   **The four that are on no list** and unchanged in v1.2 (`candidates_reads_allele.tsv`):

   | Place | HiFi reads with both anchors | ... like the assembly | Commonest version in the reads |
   |---|---|---|---|
   | chr13 paternal 28,112,400 | 17 | 1 | one C fewer in a homopolymer (10 reads) |
   | chr17 maternal 9,027,454 | 32 | 0 | one T more (22 reads) |
   | chr22 paternal 35,690,655 | 19 | 0 | 11 bases more in a GA repeat (10 reads) |
   | chr15 maternal 46,844,785 | 26 | 0 | 17 bases fewer in a GA/GAAA repeat (10 reads) |

   They are candidates, of the same kind as the listed ones, and were not checked against
   ONT or Element reads.
6. **Ambiguous segments are short.** Of the 142 segments that a depth factor explains as well as
   a whole copy, 105 are dips of one or two windows in which two-copy k-mers fall to about half
   their depth. None is both long (8 windows) and deep (1.5-fold), which is why such segments are
   reported when they occur in another assembly.

The same floor shows in the consensus quality that short-read k-mers assign to this genome:
0.018% of its tracked k-mers have no read, which reads as QV 52.4 for an assembly whose real
quality is far higher. Short-read k-mer QV saturates near there at this depth.

## Decisions

- **D5:** fixed as written in NOTEBOOK entry g; result 0.11%.
- **D10 (change of course): not needed.** After masking and classification the truth genome
  carries tens of high-confidence issues per haplotype, not thousands.
- **D1 (read source):** short reads detect; long reads adjudicate, by k-mer content, the two
  things short reads cannot decide: whether assembly sequence that no short read holds is real,
  and whether a join is supported.
- **D6:** copy-number segments that a depth factor explains as well as a whole copy
  (`CN_OR_DEPTH`) are issues only with support (k-mers of their own haplotype, the other copies
  of their sequence, another channel, or length and depth); unsupported runs are issues only with
  a missing allele from the reads or with long-read evidence.

## Files

- `analyse.py`: the tables above from a run directory.
- `issues_by_class.tsv`: issue counts by class, confidence, role and stratum.
- `issues_by_class_with_hifi.tsv`: the same for the run with `--long-read-bam`.
- `cn_candidates.tsv`: the copy-number segments labelled `CN`, with per-class copies.
- `hifi_at_short_losses.tsv`: per short copy-number loss, the haplotype's unread k-mers and how
  many HiFi reads hold them (made by `hifi_at_short_losses.py`, the hand-made first version of
  the long-read audit).
- `unsupported_runs_hifi.tsv`: every audited run with both long-read verdicts, known-issue and
  v1.2 flags (made by `annotate_runs.py`, which prints the table of item 5).
- `not_carried_reads_allele.tsv`, `candidates_reads_allele.tsv`: for the runs long reads do not
  carry, the assembly's sequence and the reads' commonest version between two anchors (made by
  `reads_allele.py`).
- `ambiguous_by_length.tsv`: ambiguous segments by length and depth, the count behind the
  long-and-deep rule (made by `analyse.py`).
- `litmus.txt`: output of `scripts/sd_boundary_litmus.py`.
