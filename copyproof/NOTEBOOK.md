# copyproof notebook

A running log of what was tried, what was seen and what changed, newest entry first. The guide proposes; this log records what the data decided.

## Entry template

```
### YYYY-MM-DD: short title
Question:  what this was meant to answer (link the GUIDE section or decision, e.g. D2)
Run:       commands or script; inputs (ids from resources.tsv); tool versions; where outputs live (outside git if large)
Seen:      the numbers that matter, with units and denominators
Decision:  what changes in the plan, or "no change"
Next:      the next step
```

## Open decisions

| # | Decision (GUIDE §9) | Status | Decided in |
|---|---|---|---|
| D1 | Read source for copy number | short reads detect; long reads adjudicate by k-mer content (`--long-read-bam`). At a false join whose other haplotype is right they adjudicate by alignment instead, from the aligner's secondary records | h, k, q |
| D2 | k-mer size | 31 throughout; other sizes not tried | a |
| D3 | Coordinate frame | the release's own sequences; nothing is projected to a reference | a |
| D4 | Segmentation | HMM over 1-kb windows with additive and depth-factor states, then classification per segment with one shared regional factor | e, i, j |
| D5 | Negative-control denominator | fixed in g before the run; result 0.11% | g, h |
| D6 | Thresholds and score | rules with log-likelihood margins; changed on rungs 1 to 4, then frozen at commit 1465496 and tested on NA19240 with the scores fixed beforehand: all seven declared comparisons hold, without a truth genome. Every value and its origin: `docs/thresholds.md` | e, h, i, j, l, m, n, p, s |
| D7 | Non-trio assemblies | S1, S2 and S4 need no parents. S3 needs both, readable; it is specific for misplaced pieces and weak for switches under 50 kb, which S2 finds without parents | c, h, l, m |
| D8 | Reuse or build | built one engine (`cpk`) | a |
| D9 | Cost | measured: about 5 minutes to count 36x of FASTQ (30 from CRAM), 5 to 10 minutes per assembly after that when the disk keeps up. A scan ends by writing a 10- to 13-GB file; on the USB disk used here that step took 1 to 20 minutes | b, l, m, o |
| D10 | Change of course | not needed: 25 high-confidence issues per haplotype on a truth genome | h, j |

## Baselines to reproduce

These values come from the pilot work in this repository. M0 and M1 should reproduce them before anything new is built on them.

| What | Value | Source |
|---|---|---|
| HG01943 read copy number in the duplicated segment (k-mers; CRAM index) | 3.13; 3.05 | report F4; `kmer_cn_track.py`, `crai_cn.py` |
| HG01943 assembly copy number in 10-kb bins | 87 bins at CN 3, 17 at CN 2, 3 at CN 4 | report F4 |
| HG01943 junction 31-mers in HPRC R2 | 0 of 28 | `regression_tests.md`, R1 |
| HG01943 piece B hapmers, maternal vs paternal | 5,638 vs 0 | R1 |
| HG01943 joins spanned by any read | 0 of 4 | R2 (T1); report F4 |
| HG01943 junction reads ending at the two ends of the maternal copy | 31 and 30 | report F4 |
| R2 targets hit: HMM-Flagger HiFi, conservative / full | 0/5 / 2/5 | R2 |
| R2 targets hit: HMM-Flagger ONT, conservative / full | 0/5 / 2/5 | R2 |
| R2 targets hit: NucFlag | 2/5 | R2 |
| NucFlag genome-wide MISJOIN calls, HG01943 paternal / maternal | 1,054 (about 361 per Gb) / 1,220 | R2 |
| HG01103 junctions in the trio assembly FASTA / seen via the published BAM | 16 of 16 / 5 of 16 | R3, R4 |
| HG01103 lost segments: paternal/maternal haplotagged depth | 0.42–0.46 (kept segments 0.95–0.97) | R5 |
| Case driver run time at 2 threads | HG01943 9 min 48 s; HG01103 8 min 2 s | `regression_tests.md` |
| Known low-confidence intervals in hg002v1.1 | 3,172 intervals, 38.8 Mb | resources.tsv: `hg002_known_issues_v1.1` |
| Samples with both parents in NYGC 30x | HPRC R2 121 of 223; HGSVC3 30 of 63 | `../data/sample_overlap.tsv` |

## Log

How this log was kept: entry g was written before its experiment, as the guide asks for a decision that must not
follow the result. Entries a to f and h to l were written in the evening of the same day, from the working record
of the session, after the experiments they describe; their numbers were not all recomputed when they were written.
Entries m, n and o were written within the hour of their experiments. Entry p, like g, was written before its
experiment; what was added to it afterwards is marked with the time. Entry s, its result, was written within the
hour of the runs. Entry t was written within the hour of its tables, which were made from finished runs.
The experiment folders are regenerated from the run directories by the scripts named in each, and are the
reference where a number here and a number there differ.

### 2026-11-14 (t): missing alleles that neither parent's reads carry, tabulated from finished runs

Note:      written within the hour, from the tables of the record
           (`experiments/2026-11-14_missing_alleles_by_parent`). Nothing was run again. No rule was made and no
           code under `src/` or `engine/` changed. Made after the held-out test and from its runs: nothing
           here is held out.
Question:  entry s saw that 922 of NA19240's 1,773 strong missing-allele calls are in neither parent's reads
           (HG002: 71 of 996) and did not act on it. Is that NA19240 alone, did the parents' reads miss
           inherited sequence, and what do the finished runs say about these calls?
Run:       `tables.py DATA OUT` of the record, 4 s. Six runs with the child's and both parents' reads: HG002's
           `hprc_year2` and the truth genome, HPRC release 2 of HG01943 and of NA19240, and the truth genome
           and NA19240's release with planted errors. The three runs of 13 November were not made again;
           `src/copyproof/missing.py`, which writes the table read here, has one commit, older than every run.
Seen:      strong missing-allele placements in neither parent's reads: 71 of 996 on `hprc_year2`, 34 of 74 on
           the truth genome, 545 of 1,701 on HG01943, 922 of 1,773 on NA19240 (38, 18, 276 and 460 events).
           The split is all or nothing: 1,561 of those 1,572 placements have no k-mer in either parent's
           database.
           The parents' reads did not miss them. Every read database keeps a k-mer from its second sighting.
           The six parents have 108.5 to 141.0 Gb of reads and the three children 107.4 to 110.6 Gb. Planted
           errors made 4,677 strong placements in the two planted runs and all are in a parent's reads (none
           in neither: at most 0.064%, one-sided 95%).
           Read counts: a median of 8 to 8.5 where the placements made by planted errors have 14; 798 of
           NA19240's 922 lie below the lower quartile of those (11). The counts begin at the threshold of 7.
           Shape: the same length as the release at all but 5 of the 1,572; a median of 7 tracked k-mers,
           which is one changed base; both strands; 36 of the 1,572 within 10 kb of another, against 1,474
           of the 2,972 placements in a parent's reads.
           Held elsewhere: by the truth genome at 51 of 71 on `hprc_year2` (in a parent's reads: 906 of 925);
           by one or both HGSVC3 assemblies at 621 of 922 on NA19240 (774 of 851). Those two assemblies,
           built from the same reads, disagree on 259 of the 922 and on 61 of the 851.
           On `hprc_year2` the read count marks the placements that the truth genome does not hold about as
           well as the parents do: of the 39, 30 have fewer than 9 reads and 20 are in neither parent's reads.
           In the issue lists: 906 of NA19240's 1,258 high-confidence `MISSING_ALLELE` issues, 539 of
           HG01943's 1,085, 71 of 734 on `hprc_year2`, 34 of 62 on the truth genome.
Decision:  none on the method. The tables fit sequence changes that are in part of the child's cells and in
           neither parent, as changes arising in a cell line would be; no table shows that they are that.
           Two sentences were made exact: the README's 96% for strong missing alleles now names HG002, and
           `docs/method.md` says that "few in neither parent" counts one-sided events.
Next:      for a review to decide: whether such calls get a wording or a role of their own when parents are
           given, and whether a rule should use the parents, the read count or both. HG002, HG01943 and
           NA19240 have all been looked at for this, so a rule needs another person to test it. What the
           calls are could be read from the changed bases, from the long reads each release was built from,
           or from a second library; none was looked at.

### 2026-11-14 (s): rung 6, the held-out test on NA19240

Note:      written within the hour of the runs, from the tables of the record
           (`experiments/2026-11-14_na19240_heldout`). Entry p fixed the test beforehand; nothing in it was
           changed afterwards except the additions it marks with their times.
Question:  with every threshold fixed beforehand, do the results hold on a second person (GUIDE §8; D6)?
Run:       `scripts/run_na19240_heldout.sh`, and `plant.sh`, `alt_settings.py`, `tables.py` and `after.py` of
           the record, with `src/` and `engine/` at commit 1465496. Inputs of
           `scripts/manifests/05_na19240_heldout.tsv`, each assembly and CRAM checked against its published MD5;
           one file of the HGSVC3 hifiasm assembly had arrived with 198 wrong bytes and was repaired. The
           pre-registration was committed (a82fdbe) and pushed at 10:29; the scored run started at 10:30:45 and
           took 1,167 s. The run with long reads, the planted genome, the Verkko assembly and the alternative
           settings followed one at a time until 11:41. Outputs: DATA/work/rung6 (95 GB); k-mer databases under
           DATA/kmerdb (72 GB in 11 files).
Seen:      all seven declared comparisons hold (`verdicts.tsv`). Score 1: 46 of 81 high-confidence copy-number
           issues confirmed by both other assemblies, predicted 74 of 119 (p = 0.27); by the Verkko assembly
           alone 75, by the hifiasm assembly alone 53, by neither 6. Score 2: 64 of 76 surplus copies of 5 kb
           and more found, predicted 117 of 127 (p = 0.066); 2 of 14 collapses, predicted 1 of 3. Score 3: 47
           of 48, 28 of 31, 41 of 43 and 8 of 12, none below HG002's share.
           What "holds" is worth. Against the second HG002 pair, 97 of 119, score 1 is lower (p = 0.0001).
           The difference is the second judge. The HGSVC3 hifiasm assembly holds 7.17 Gb; of the 2,813
           stretches of 5 kb and more where it alone holds more copies than the release and no issue lies, the
           reads side with the release at 2,777, 1.38 Gb. It confirms 44 of 69 surplus copies where the Verkko
           assembly confirms 64.
           Of the 12 surplus copies not found, the classifier saw 7 and did not call them (3 judged a depth
           change, 4 ambiguous). Of the 12 collapses not found, 11 have no segment, and at 10 the reads stand
           at the release's count: this proxy cannot judge collapses.
           Missing alleles: one or both assemblies hold 1,395 of 1,773 strong calls (79%; HG002 92%). Calls
           whose k-mers a parent's reads carry are held at 91% (HG002 93%). 922 calls are in neither parent's
           reads (HG002 71), with a median of 8 reads where one copy gives 13.6.
           With long reads: `scripts/eval_assemblies.py` scores issues of channel S1 only and so leaves out
           the 38 that the audit marked S1+S4. Counted in full: 79 high-confidence copy-number issues, 44
           confirmed by both and 72 by the Verkko assembly. The audit made 3 of the 5 high-confidence surplus
           copies it looked at low-confidence; the Verkko assembly confirms all 3.
           The two candidate values: with `min_seg_class` 10 and `margin` 2, 84 high-confidence issues where
           the shipped values give 78, 79 confirmed by the Verkko assembly where they give 72, and 67 of 76
           surplus copies found where they give 64.
           The added run: on the Strand-seq-phased Verkko assembly the hapmer audit calls 65 switch segments,
           13 of them with 10 sites or more. HPRC's trio-phased haplotypes confirm 7 of the 13 and contradict
           none; of the 52 weaker ones they confirm 16 and contradict 8. Two stretches of 12.4 Mb, each at the
           end of a sequence, make up 24.8 of the 27.2 Mb.
           Costs: the three read counts took 2,721 s, 2,450 s and 758 s for 109, 141 and 109 Gb. The first two
           ran while a CRAM was being downloaded to the same drive; which load cost the time was not measured.
Decision:  D6: the thresholds stay as shipped. Frozen, they held on a second person by every declared
           comparison. The test has no truth genome: its two measures of the share of true issues give 57% and
           93%, as they gave 62% and 85% on HG002 where the truth gave 96%. The two candidate values looked no
           worse on the new sample; changing them is a separate decision and was not made here. No code under
           `src/` or `engine/` changed.
           Seen and not acted on: the scoring script's channel filter on runs with long reads; the long-read
           audit on an assembly made from the same reads; missing-allele calls in neither parent's reads.
Next:      not decided here. The scan of whole releases, the guide's rung 6, was not started. A second person
           with a truth genome would settle what this test could not; the candidate found, I002C, has MGI
           short reads and no parental short reads (`resources.tsv`).

### 2026-11-14 (r): the 69 unconfirmed issues of rung 4, placed in the truth genome and sorted

Note:      written within the hour of the runs, from the tables of the record. Nothing under `src/` or `engine/`
           is changed (the code is frozen at commit 1465496 for the held-out test). No threshold was changed and
           no call was relabelled: the scoreboard of rung 4 stands as it was.
Question:  are the high-confidence copy-number issues that the truth genome does not confirm the tool's errors,
           errors of the truth genome, or differences between the DNA that was read and the genome that was
           assembled (review of 2026-11-14, point 3.2)?
Run:       in `experiments/2026-11-13_hg002_assemblies/`: `place_unconfirmed.py` (`cpk extract` on each of the 13
           assemblies; `cpk locate --query` of the 69 calls, 8.6 Mb, in T2T-HG002 Q100 v1.1, 1 k-mer in 4, up to
           64 occurrences; each call's k-mers from its run's track), `sort_unconfirmed.py`, then `locus_depth.py`
           on six loci (48 stretches: the track of the trio run of entry m, and the public HiFi alignment read by
           region). Record: the section "The 69 unconfirmed issues, sorted" of that folder. Outputs:
           DATA/work/hg002_assemblies/unconfirmed (0.6 GB). The steps that load a genome ran for 25 minutes,
           8 of them in one stall of the drive, and waited 62 minutes behind the held-out test's jobs.
Seen:      All 69 are placed (1,283 placements). Where the benchmark's own alignment covers a call (52), it
           overlaps a principal placement of 50.
           27 of 69 are DNA of which a cell holds another number than a genome lists. 15 are Epstein-Barr virus
           contigs in five assemblies: the truth genome v1.1 holds none of their k-mers, and the reads hold 31 to
           37 copies. 2 are a contig that holds the mitochondrial genome twice (1,113 copies in the reads). 10
           are the pair of calls 57 to 64 kb into chromosome 1 that entries h and o left unexplained: the stretch
           is an insertion of mitochondrial sequence, 56% to 60% of its k-mers also occur in chrM and are read
           at about 1,100 copies, and its own k-mers are read at 1.07 of the truth genome's count in short reads
           and 1.03 in long reads. "Four more copies" was the cap of the window statistics, not a count.
           24 of 69 are real differences from the truth genome: it holds another number of copies of at least
           half of the call's k-mers for 22 (a fifth to a half for 2). The scoring rule, which wants such
           k-mers to be the majority in windows covering half of a call, fell short for three reasons: 13 calls
           lie in sequence repeated more than 8 times, 4 are mostly gap (63% to 86% N), and 5 are wider than
           the difference (it covers 40% to 48% of them).
           18 of 69 have only the reads behind them; all are expansions; 15 lie on acrocentric short arms. 12 are
           two loci on the paternal chromosome 13: an HSat1B block at 7.26 Mb that begins a 345-kb duplication
           (10 calls in ten assemblies) and an 11-kb four-copy stretch at 10.04 Mb (2 calls, and the truth
           genome's own run).
           The test that told the cases apart: count the locus's k-mers in the short reads of HG002 and of both
           parents and in HG002's long reads. At the HSat1B block the long reads hold the truth genome's count
           (0.91 to 1.09 of it) and the short reads of all three people fall short (0.37 to 0.81): short reads
           of this kind under-read the satellite, and the classifier took a shortfall that differs by class for
           a copy. The
           same at chr9 42.4 Mb (60 k-mers at 0.42 to 0.48 in three people, 0.98 in long reads) and, mildly, at
           chr7 58.7 Mb. At chr13 10.04 Mb and at chr22 3.9 Mb the long reads fall short too (0.83 to 0.99; 0.73
           to 0.93). Every copy there has long reads running through it, and over all 29 short-arm stretches the
           long reads' alignment depth is a median 0.89 of that over ordinary sequence.
           None of the 52 nuclear calls has half of a placement on a known issue (the most is 27%); v1.2 changed
           four positions inside placements, by 1 to 14 bases. 6 of 27 loci recur in two or more independent
           assemblies and hold 42 of the 69 calls.
           Three things went wrong on the way. The wait for the held-out test's jobs matched its own command
           line and waited on itself for a minute. A second placement pass for high-copy k-mers wrote 3.7 GB
           and was not needed, because the track already holds each k-mer's count in the truth genome; it is
           removed. And the chromosome 13 calls were first read as one copy missing from the short reads' DNA;
           the parents' reads and the long reads showed otherwise.
Decision:  no change to the method while the code is frozen. Nothing here is an error of the truth genome that
           could be shown; two loci are written up for its curators as open questions, with what is missing.
           For after the held-out test, as suggestions and not as changes: sequences whose k-mers are read at
           tens to thousands of copies (virus, mitochondrion, and nuclear copies of mitochondrial sequence)
           should be named as such and not called collapsed; a count at the cap should be reported as "at
           least"; and with `--long-read-bam`, the long reads' count of a loss's k-mers could lower its
           confidence when it equals the assembly's, which is now done only where no short read holds the
           k-mers at all.
Next:      the two loci need reads of a third kind (ONT, Element) counted the same way. A score that also counts
           k-mers above multiplicity 8 would say how many of the 24 real differences rung 4 undercounts; not
           done.
### 2026-11-14 (q): rung 2 with long reads: what the audit says at planted false joins

Note:      written within the hour of the runs, from the tables of the record. The code under `src/` is
           frozen at commit 1465496 for the held-out test (rung 6) and is unchanged here.
Question:  given real HiFi reads, what does the long-read audit say at the false joins planted in the truth
           genome (entry i, "not done"; review of 2026-11-14, point 3.4; D1)?
Run:       `scripts/run_planted_hifi.sh`. The reads that the public HiFi alignment places within 50 kb of an
           edit or of its homologous stretch (222 stretches, 33.8 Mb; 80,615 reads, 1.30 Gb) were fetched by
           region and aligned again to the planted genome with minimap2 2.30 (`-ax map-hifi`). The pipeline ran
           with `--long-read-bam` on links into DATA/work/plant/run35, once with that alignment and once without
           its secondary records. Scoring: `experiments/2026-11-13_planted_errors/long_read_audit.py`. Record:
           the part "With long reads" of that folder. Outputs: DATA/work/plant/hifi (2 GB), run35_hifi and
           run35_hifi_no_secondary. The fetch and the alignment were first run by hand, with the commands that
           are now in the driver; the driver ran both pipeline runs and the scoring. Fetch 7 minutes, alignment
           4, each pipeline run 3.
Seen:      `FALSE_JOIN`, all of high confidence, at 23 of 24 inversion breakpoints and at 21 of 36 joins left by
           deletions: 16 of 20 from 10 kb, 5 of 8 at 5 kb, 0 of 8 at 2 kb. All 49 `FALSE_JOIN` issues lie at a
           planted adjacency (5 are junctions of tandem duplications), all 86 issues the audit made or changed
           lie within 5 kb of an edit, and the scoring of rung 2 changes in the nine rows of the deletion's join
           and in no other. The 8 runs of N of the inversions rejoined through N have no linking read: 8
           `UNSUPPORTED_JOIN`, of low confidence, because nothing else is wrong within 20 kb of them.
           The misses have two causes. The audit asks only about runs of unread k-mers whose flanks short
           reads cover: that leaves out 1 breakpoint, 7 deletion joins and 32 of 39 duplication junctions. And
           at the 2-kb deletions 405 alignments run through the 8 joins, each with the deleted 2 kb as an
           insertion, 402 of them secondary records and none with a mapping quality of 20 or more; the rule
           counts them as reads through the join (2 `BASE_ERROR`, 3 not reported among the 5 asked about).
           The calls rest on alignment, not on content. No long read holds an unread k-mer at any of the 59
           joins asked about, but "not carried" needs 5 reads with their sequence there and the median is 2,
           because a read that crosses a false join gets its primary record on the other haplotype. Content is
           "unclear" at 51 of 59, and the rule "no alignment through, 3 or more clipped" decides. Of 582
           records clipped at the inversion breakpoints, 534 are secondary, all of reads whose primary record
           is on the other haplotype. Without the secondary records: `FALSE_JOIN` at 8 of 24 breakpoints and 3
           of 36 deletion joins.
           The aligners differ on the same reads: minimap2 wrote 85,778 secondary records where the public
           Winnowmap alignment has 21,697, and put 9,342 primary records on another sequence.
           Away from the edits the audit asked about 4 runs in the fetched stretches; long reads hold their
           k-mers, as in the truth genome's run. Four is too few to call a false-flag rate.
           Beside the experiment: a run's `run.json` does not name the long-read alignment, because the
           pipeline does not record `--long-read-bam`; the address of the public file is in `resources.tsv`
           (`hg002_q100_hifi_bam`). The `run.json` of these two runs does name their alignment: they were made
           with one more line in `pipeline.py` that writes it, and that line is not kept, because the code is
           frozen. It changes nothing else that a run writes. The drive paused three times, for 10 to 15
           minutes each, while other jobs wrote to it (entry o); the steps waited and finished.
Decision:  no threshold or rule changed. D1 is qualified: at a false join whose other haplotype is right, long
           reads adjudicate by alignment, from the aligner's secondary records, and not by k-mer content. The
           documents now say that the alignment file must keep its secondary records.
Next:      not done, in the order they would matter. (1) The same test on a Winnowmap alignment and with ONT
           reads. (2) After rung 6, when the code may change again: write `--long-read-bam` into `run.json`,
           and three changes to the audit that these numbers suggest and that were not tried: take the sequence
           of a read from its primary record wherever that lies, so that content can speak at a false join; do
           not count an alignment as running through a run when it holds a long insertion or deletion there,
           or when its mapping quality is under 20; let the two unlinked joins of one inversion vouch for each
           other. Each of the three needs the truth genome's count of false flags again. (3) A false-flag
           count for the minimap2 alignment on unedited sequence: fetch the reads at the 1,123 runs of rung 3
           and compare the verdicts with those of the public alignment.

### 2026-11-14 (p): rung 6 fixed before the run

Question:  with every threshold fixed beforehand, do the results hold on a second person (GUIDE §8 rung 6; D6)?
Run:       none yet. Written at 07:34 on 2026-11-14, when the downloads had started and before any read of the new
           sample was counted or any program had opened its assemblies. The sample was chosen from two
           candidates; the other was I002C, a second genome with a near-truth assembly, whose public short
           reads are MGI and whose parents' short reads are not public.
Seen:      n/a.
Decision:  D6.
           Sample. NA19240, a trio child of HPRC release 2 on which copyproof has not been run. Release under
           test: `NA19240_pat_hprc_r2_v1.0.1.fa.gz` and `NA19240_mat_hprc_r2_v1.0.1.fa.gz` (hifiasm 0.19.9, trio).
           Reads: the NYGC 30x CRAMs of NA19240, her mother NA19238 and her father NA19239, each counted once.
           Long reads: HPRC's HiFi-to-assembly alignment, read remotely. Manifest:
           `scripts/manifests/05_na19240_heldout.tsv`.
           What is frozen. The code of commit 1465496 with every default it carries (`docs/thresholds.md`). The
           drivers refuse to run when `src/`, `engine/src` or `scripts/run_sample.py` differ from it. Of the two
           values that were not the best on the old data, the shipped ones are tested (`min_seg_class` 20,
           `margin` 3), because every reported number uses them.
           Evidence. No truth genome exists for this person. Two other assemblies of her do: the Verkko and the
           hifiasm assembly of HGSVC3, made at another centre from other reads. Their exact k-mer counts label
           every tracked k-mer of the release as a truth genome's do in rung 4, with one rule added, called
           `both`: a k-mer is contradicted when the two assemblies hold the same number of copies of it and
           that number differs from the release's. All else is rung 4's. Only k-mers that the release holds at
           most 8 times count. A discordance is 2 or more 1-kb windows (joined across one window), each with at
           least 30 such k-mers of which half or more are contradicted in one direction. An issue is confirmed
           when discordances of its direction cover at least half of it. A discordance is found when
           copy-number issues of its direction cover at least half of it, and is assessable when under half of
           it lies in the run's mask. Script: `scripts/eval_assemblies.py`, tested on the synthetic trio
           (`tests/test_eval_assemblies.py`).
           This is a proxy, not truth. Two assemblies can share an error, and the hifiasm one more readily,
           because the release is a hifiasm assembly too. An issue that neither confirms is not thereby false.
           So every count is also given for the Verkko assembly alone, the hifiasm assembly alone, and either.
           The scored run. The child's reads and both parents, no long reads (`short`), because neither the
           calibration below nor rung 4 has long reads. The run with long reads (`long`) is reported beside it.
           Driver: `scripts/run_na19240_heldout.sh`.
           Score 1. Denominator: every high-confidence copy-number issue of the `short` run (channel S1; class
           COLLAPSE, EXPANSION or UNSUPPORTED; any role). Count: those confirmed by `both`.
           Score 2. Denominator: the assessable discordances by `both` of 5 kb and more, separately for the two
           directions. Count: those found.
           Score 3. Errors planted in the release with the design of rung 2: the same kinds, sizes and number
           per cell, seed 20261113, 3 Mb apart. Placed on sequences of at least 12 Mb (the release has no
           chromosome names; the sample is female, so X is not set apart), outside HPRC's satellite annotation,
           outside HPRC's own Flagger and NucFlag calls (in place of a list of known issues) and outside the
           copy-number mask of the `short` run. Duplications come from HPRC's annotation, the homologous
           stretch from an alignment of the two haplotypes (minimap2, asm5). Run with the child's reads and
           both parents, scored by `scripts/eval_planted.py`, unchanged, against the `short` run. Counts, each
           from 5 kb: false duplications (kinds dup and extra, signal EXPANSION); deletions as a collapse on
           the homologue (del, COLLAPSE); deleted or overwritten stretches as a missing haplotype (del and
           swap, MISSING_HAPLOTYPE); and swaps as a phase switch (swap, PHASE_SWITCH, every size).
           Prediction. For scores 1 and 2, the same rule applied to HG002, where the truth is known: the HPRC
           year-2 assembly (`hprc_year2`), run with HG002's reads and parents, judged by `v2_trio` (Verkko) and
           `hprc_year1` (hifiasm, other reads), with the truth genome's verdict beside it
           (`experiments/2026-11-14_na19240_heldout/calibrate_hg002.sh`). It runs after this entry and its
           counts are added below before the first run of NA19240. A second pair (`verkko_asm_thic`,
           `hprc_year2_v2_thic`) shows how far the proxy moves with the assemblies chosen and predicts nothing.
           For score 3, the counts of HG002's planted genome with its parents: 46 of 47, 25 of 28, 34 of 40
           and 7 of 12.
           Verdict. A score holds when NA19240's share is not below the predicted one by a one-sided Fisher
           exact test at 5%. It is reported with its counts either way. Seven comparisons are made (score 1;
           score 2 in two directions; score 3 in four rows), with no correction for their number. HG002's other
           assemblies are not HGSVC3's, so the prediction is an anchor, not a guarantee.
           Reported without a prediction. Issues by class and confidence beside HG01943's. Strong missing-allele
           calls whose k-mers the other assemblies hold. Assembly k-mers that no read holds and the other
           assemblies lack. Misplaced pieces, phase switches and joins. Overlap of the issues with HPRC's
           Flagger and NucFlag calls. Scores 1 and 2 for the `long` run.
           The two candidate values. The segments of the `short` run are classified again with `min_seg_class`
           10, with `margin` 2 and with both, as in the sensitivity record (without the lookup of other copies
           and without support from other channels), and scored by `both`. This is a second, declared result.
           It changes neither the scored run nor the shipped values.
           Breaking the freeze. If code under `src/` or `engine/` has to change after the first output of
           NA19240 is seen, that is a new entry and the test is no longer held out; the record will say so.
           Not committed beforehand: commits to this branch are made on request only. The order rests on the times in
           this entry, on the logs under DATA/logs and DATA/work/rung6, and on the checksums of the scoring
           script and drivers in `experiments/2026-11-14_na19240_heldout/prereg.sha256`.
           Prediction, added at 07:53 on 2026-11-14, after the calibration and before any program had read a
           file of NA19240 other than to download it or check its checksum (`predictions.tsv` and the `hg002_*`
           tables of the record). By the truth genome the calibration run gives rung 4's counts again: 114 of
           119 high-confidence copy-number issues confirmed, 151 of 233 assessable discordances found. By the
           proxy, on the same 119 issues: `both` 74, the Verkko assembly alone 101, the hifiasm assembly alone
           82, either 103. `both` confirms no issue that the truth does not; the Verkko assembly alone confirms
           2 of the 5 that the truth does not. Collapses are where the proxy fails: of 24 high-confidence
           collapse issues the truth confirms 22, the Verkko assembly 9, the hifiasm assembly 2 and `both` none.
           So score 1 will speak for surplus copies almost alone.
           Score 1: 74 of 119. Score 2: 1 of 3 where the other assemblies hold more, and 117 of 127 where they
           hold fewer. With three discordances the first of these can show next to nothing; it stays because it
           was declared. Score 3: 46 of 47, 25 of 28, 34 of 40, 7 of 12.
           The second pair gives `both` 97 of 119. The proxy thus moves by 19 points with the assemblies chosen,
           and the declared prediction is the lower of the two anchors: "holds" will say little by itself, and
           the record will set NA19240's share against both. Strong missing-allele calls on HG002, for the same
           comparison later: the truth holds the k-mers of 957 of 996, both other assemblies of 706, either of
           913.
           Two points from a review of the pull request, added at 08:00 on 2026-11-14, while the first CRAM of
           NA19240 was still downloading and none of her reads had been counted.
           The hapmer channel. The release was phased with the parents' k-mers (hifiasm, trio), so an audit with
           the same parents is not an independent test of its phasing. It can show what went wrong after
           phasing, such as a piece scaffolded into the other haplotype (HG01943), and nothing more: few
           switches there will not mean that the phasing is right. Both HGSVC3 assemblies were phased with
           Strand-seq and without the parents (the README files of their two batches say so), so the same audit
           on one of them would be independent. That run is not part of this test unless it is added.
           Counts, not only shares. Score 1 is reported as the count confirmed and the count not confirmed, with
           the share, and the issues not confirmed are listed one by one. Rung 4's share of 99% rests on one
           assembly with 3,390 of the 4,648 calls; without it the share is 95% (1,199 of 1,258).
           A run added on request, written down at 08:05 on 2026-11-14, still before any read of NA19240 was
           counted: the hapmer audit where it is independent. The HGSVC3 Verkko assembly, every file of it, is
           run with the same reads and parents at the frozen thresholds, without long reads and without the
           missing-sequence channel. It was phased with Strand-seq, so the audit here tests phasing that the
           parents did not help to make. Each call is then held against HPRC's trio-phased haplotypes, whose
           exact k-mers ride along as two label sets. Among the tracked k-mers in a call that only one HPRC
           haplotype holds, those of the parent the call names must number at least 20 and at least three
           times those of the sequence's own parent: then the call is confirmed. The reverse makes it
           contradicted, anything else unclear. The other way round, the assessed sequences are cut into 50-kb
           bins. A bin follows the other parent by the same rule of 20 and three to one, and is flagged when
           calls cover at least half of it; a bin thus speaks only for switches that fill most of 50 kb.
           Reported: the sequences assessed and whether hapmers and HPRC's k-mers name the same parent for
           each; switch segments and misplaced pieces confirmed, contradicted and unclear; the bins that follow
           the other parent and how many of them are flagged. No prediction: nothing like it was measured
           before (on HG002's truth genome, where every switch is false, the audit left 24 segments). HPRC's
           haplotypes were phased with these parents' k-mers from other reads: another assembly, not another
           pedigree. Only the hapmer channel of this run is part of the record. Script:
           `experiments/2026-11-14_na19240_heldout/hapmers_against_hprc.py`, tried on the synthetic trio.
           Committed beforehand after all, added at 10:31 on 2026-11-14, while the scored run was in its first
           minutes and before any output of it had been read. On request, this entry and the 28 other files of
           the pre-registration were committed as a82fdbe at 10:27:36, merged with the branch as 8f5b972 and
           pushed. GitHub's activity record gives the push 16:29:27 UTC, which is 10:29:27 here. The scored
           run was held for it and started at 10:30:45. The files of a82fdbe are those of `prereg.sha256`,
           so the order now rests on GitHub's time for the push as well as on the local times above.
Next:      NA19240.

### 2026-11-13 (o): six more assemblies; the disk measured while a scan waited on it

Question:  does rung 4 hold on more assemblies, among them drafts of the truth genome, where almost nothing is
           wrong? And a second question about the machine: is the data disk slow now?
Run:       `scripts/run_hg002_assemblies.sh` for `hg002v0.9`, `hg002v1.0.1`, `hprc_curated`, `verkko_asm_thic`,
           `v2_trio`, `hifi_q28_trio_hic` and `lc24_nopolish`, one at a time; `scripts/scoreboard.py`, which
           now also writes out the sums the record quotes;
           `experiments/2026-11-13_hg002_assemblies/unconfirmed.py`. For the disk: cores used by the running
           `cpk scan`, growth of its track file and MB moved by the device, every 5 s; a stack sample of the
           process; write tests of random data flushed with F_FULLFSYNC.
Seen:      Assemblies, thirteen now, which is every one that was downloaded: 4,579 of 4,648
           high-confidence copy-number issues confirmed (99%; 1,199 of 1,258 without `verkko_2023`). Surplus
           copies of 5 kb and more: 73% to 79% found. Collapses: 25% to 54%. Under 5 kb: 9% to 18%. The
           unconfirmed ones are 1 to 10 per assembly whether the confirmed ones are 0 or 3,380. The three
           drafts have 3 unconfirmed each beside 4, 3 and 0 confirmed; two of the three
           are the same in every draft: the pair of calls near the start of chromosome 1 that the truth genome's
           own run makes. QV by read k-mers saturates: 51.2, 51.8 and 52.0 for drafts that the truth's k-mers put
           at 56.9, 58.7 and 59.3. `hprc_curated`, which has accession names and no haplotype labels, scores
           like `hprc_year1`.
           One more run, `hifi_q28_trio_hic`, gave the numbers of `verkko_asm_thic` to the last digit. The
           two downloads are one file (same size, same MD5): the benchmark folder holds that assembly under
           two names. It is counted once, and it is an unplanned repeat: the two runs, 43 minutes apart, wrote
           identical issue lists, segments, depth anomalies, models and scores.
           The disk: entry m says the evening's slow scans were not the disk. That was wrong, and it was a
           guess: short write tests had been fast, and memory bandwidth was blamed without being measured.
           Measured now. The track-writing step of six scans that evening took 95 s to 1,218 s for
           9 to 12 GB (95 MB/s down to 10). At 22:06 a 4-GB test wrote 690 MB/s for 2.3 GB, then 270 to 300. At
           22:10 the next scan began its track. For six and a half minutes it used 0.0 cores, all 1,766 stack
           samples of its main thread were inside `write`, the file grew by about 2 MB/s and the device moved 1
           to 3 MB/s; a new file of random data took 32 MB in 43 s. The same stalls are in the earlier logs:
           200 s to load a FASTA that otherwise loads in 28 to 48 s, and 288 s to write a table of hits.
           Then a second wrong turn. That scan was paused for five minutes, and when it continued it wrote
           the other 10.8 GB in 55 s (170 to 235 MB/s, 3.3 to 3.9 cores busy). The pause was taken for the
           cure and a watchdog was set to repeat it. On the next scan the track stalled after 5.7 GB; the
           watchdog paused it twice for five minutes, and both times the drive was still at 1 to 2 MB/s
           afterwards. Left alone, the stall ended by itself, 15 minutes after it began, and the other 6 GB
           went out at 215 MB/s. The pause did nothing that waiting would not have done.
           What the numbers support: this drive takes writes at 200 to 700 MB/s for a while, then nearly stops
           for 10 to 20 minutes, reads included, and resumes by itself. Between the end of one stall and the
           start of the next, 17 to 18 GB were written. A scan writes 12 GB, so nearly every run meets one
           stall, and an assembly with its scoring took 16 to 26 minutes, most of it the stall. A short
           write test shows none of this, and what the drive does during a stall was not established.
Decision:  no change to the method. The costs in the README hold when the disk keeps up, and it now says that
           a scan ends by writing a 10- to 13-GB file in one stream. Entry m and D9 are corrected.
Next:      none for rung 4 on this sample. Rung 5 is other samples.

### 2026-11-13 (n): a review pass over the day's claims

Question:  can a reviewer follow the findings and the reasons for the decisions from what is in the repository?
Run:       every headline number was traced to a table and the table to a script; see `experiments/README.md`.
Seen:      Six gaps. (1) Several numbers came from analyses typed once and not saved. They are scripts now and
           were run again: the known-issue and v1.2 comparison (`annotate_runs.py`), the planted breakpoints (in
           `eval_planted.py`), the counts behind the long-and-deep and short-call rules (`analyse.py`), the missed
           collapses and cluster precision (`details.py`), the HG01943 tables (`tables.py`), the jellyfish
           comparison (`compare_jellyfish.py`). All reproduced; two corrections: 46 events in neither parent, not
           "about 45", and the jellyfish table moves by a few k-mers between runs.
           (2) The classifier's history was prose. The classifier before its revision is kept in `scripts/dev/`
           and both are run on the same cached segments (`experiments/2026-11-13_threshold_sensitivity/`).
           (3) Entries h, i and l quoted counts from before later changes; each now says so.
           (4) Most thresholds had no stated origin. `docs/thresholds.md` lists every one as fitted, set from a
           count, or chosen. The sensitivity check found most on a plateau, and two chosen values that are not
           the best on these data (`min_seg_class` 20 against 10; `margin` 3 against 2). They were left alone:
           another change on the same data would be one more, and they are the first candidates for rung 6.
           (5) Findings were coordinates without evidence. For the 119 runs of Q100 that no read holds, the
           reads' own version was extracted (`reads_allele.py`): in 110 runs reads can be compared, at most one
           read agrees with the assembly in each, and 106 differ in length (28 in a homopolymer). The four on no
           list are a homopolymer, a 1-base insertion and two GA-repeat lengths. For HG01943 the "seven more
           joins like the known four" did not hold up: at the known joins the surplus copy repeats sequence from
           beside another join of the same scaffold; at none of the seven does it, and five are at gaps of 19 to
           795 kb. The record says so now.
           (6) There was no map from claim to table to script, and this log did not say how it was kept.
Decision:  no threshold changed in this pass. The records are regenerated from the run directories and are the
           reference where they and this log differ.
Next:      rung 6 with `min_seg_class` and `margin` decided beforehand.

### 2026-11-13 (m): the hapmer channel on the truth genome; the disk was not slow

Note:      the sentence below that later slowdowns "were not the disk either" is wrong, and so is the time
           given for one scan; entry o has the measurement. The drive was fast in the short test at 20:30 and
           stalled under the scans' 10-GB writes.
Question:  what does S3 flag on a genome whose phasing is right, and does it find planted haplotype swaps (D7)?
           Taken up after a doubt that the disk was still slow.
Run:       `scripts/run_hg002_trio.sh` (parents HG004 and HG003, counted in 5 minutes each), then
           `experiments/2026-11-13_hg002_trio/tables.py`. Outputs: DATA/work/hg002_v1.1/run35_trio and
           DATA/work/plant/run35_trio.
Seen:      The disk first: it had been reported slow from a test at 18:47 (2 MB/s), and that was repeated an
           hour and a half later without testing again. At 20:30 it wrote 545 MB/s. The cancelled runs were
           restarted.
           Later slowdowns that evening were not the disk either (still 400 to 650 MB/s): three scans at once
           compete for memory bandwidth and each ran ten times slower. One scan at a time takes 2 to 5 minutes.
           S3 on the truth genome: all 46 chromosomes on the right parent, no misplaced piece among 18, but 143
           phase-switch segments (28 with 10 sites or more, 10.9 Mb) where there should be none. The largest are
           the immunoglobulin loci on the maternal haplotype (up to 151 sites against 7): the mother's sample is a
           B-cell line that has deleted part of the haplotype she passed on. Of the 143, 111 lie where a parent's
           depth on two-copy k-mers is outside 0.7 to 1.5 of its usual depth. Seven are a handful of sites (5 to
           14) spread over 0.7 to 1.4 Mb each.
           Planted swaps: S3 finds 1 of 4 at 10 kb, 2 of 4 at 50 kb, 4 of 4 at 200 kb; the missing-haplotype call
           of S2 finds 11 of the 12 without parents.
           S2 with parents: a rule written without data made strong one-sided missing sequence
           high-confidence when a parent's reads held the k-mers. On the truth genome a parent holds the k-mers of
           208 of 211 strong clusters. Missing k-mers in neither parent: 5 events.
Decision:  D7. Hapmer sites count only where both parents are readable, and a switch needs sites at a tenth of
           the usual density (0.1 to 0.2 give the same result on the truth genome; 0.2 loses one planted swap):
           143 segments become 24 (0.4 Mb; 8 with 10 sites or more, 4 of them still the heavy-chain locus), and
           the planted swaps are found as before. The inheritance rule for missing sequence is removed.
           HG01943 with these changes: phase switches 156 to 32, the misplaced piece still called (148 sites
           against 0), the locus score unchanged at 5 of 5.
Next:      a phase switch at an immunoglobulin locus of a cell-line trio could be labelled as such, if the loci
           were given.

### 2026-11-13 (l): rung 1, HG01943 in HPRC release 2

Note:      entry m changed the hapmer channel and one missing-sequence rule after this was written. The final
           counts are in `experiments/2026-11-13_hg01943_trio/` (2,737 issue rows, 32 phase-switch segments).
           The "7 more joins" below turned out weaker candidates than this entry suggests: see that record.
Question:  on the published assembly with public reads, are the five defects of the HG01943 locus flagged without
           flooding it (GUIDE §8 rung 1; regression test R2), and what is reported genome-wide?
Run:       `scripts/run_hg01943.sh`: NYGC 30x CRAMs of the trio, HPRC's HiFi-to-assembly BAM read remotely.
           Record: `experiments/2026-11-13_hg01943_trio/`. Outputs: DATA/work/hprc_r2/HG01943.
Seen:      First run: 3 of 5 targets (the four joins, the misplaced piece, both collapsed stretches; one false
           duplication; no junction), 1 off-target call. The junction's reads had been placed within 1 bp, but as
           a weak one-sided cluster of 4 reads, because its far side is a high-copy repeat where reads are not
           recruited; the second false duplication was an ambiguous "3 of 4 copies" right at an unlinked join.
           After two general rules (weak missing sequence is reported when its k-mers are in neither parent;
           an ambiguous copy-number segment within 5 kb of an unlinked join is reported): 5 of 5, 2 off-target,
           pass. The published QC of this assembly scores 0 to 2 of 5.
           Genome-wide: 2,867 issue rows. One misplaced piece in the whole genome, the known one (281 maternal
           hapmer sites against 0 in a paternal scaffold). 34 of 42 scaffold joins have no linking long read, so
           an unlinked join alone is uninformative; 11 lie next to another issue, the 4 known and 7 more with a
           copy too many beside them. 156 phase-switch segments (19 high-confidence), 12.7 Mb. Of 2,184 runs of
           assembly k-mers that no short read holds, HiFi reads carry 1,072 and not 659.
Decision:  an unlinked join is high-confidence only next to another issue. The two rules above. R2 passes, but it is
           not an independent test of those two rules. Their cost in this genome is 2 segments and 46 events; the
           second rule's cost on a truth genome is not measured, because that needs HG002's parents.
Next:      the 7 other joins with a neighbouring expansion are candidates of the same error class; not examined.
           HG002's parents (HG003, HG004) are downloaded and `scripts/run_hg002_trio.sh` is ready: it gives the
           false-flag rate of the hapmer channel on the truth genome and its sensitivity on the planted swaps.
           It was started and stopped: after about 450 GB of writes in a day the USB SSD dropped to 2 MB/s.
           (Run two hours later: entry m.)

### 2026-11-13 (k): long reads settle what short reads cannot

Question:  when no short read holds an assembly k-mer, is the assembly wrong or did short reads fail (D1)?
Run:       the truth genome again with `--long-read-bam` pointing at the public HiFi-to-Q100 alignment (read by
           region, nothing downloaded but slices). Outputs: DATA/work/hg002_v1.1/run35_hifi; tables
           `unsupported_runs_hifi.tsv` and `hifi_at_short_losses.tsv` in
           `experiments/2026-11-13_hg002_q100_false_flags/`.
Seen:      First by hand at the 13 stretches of 1 to 3 kb where one haplotype's own k-mers have no short reads:
           HiFi reads carry all of those k-mers at 10 (the sequence is real), few or none at 2 (chr5 paternal
           188.49 Mb, chr18 maternal 23.78 Mb, both known Q100 issues), and one has too few reads to say.
           Then for all 1,123 runs of unsupported assembly k-mers with covered flanks: HiFi reads carry the k-mers
           of 957, do not carry those of 119, and 47 are unclear. Of the 119, 115 lie on a known Q100 issue
           interval and 27 were changed in v1.2; of the 957, 74 and 22. Four not-carried runs are on no list
           (chr13 paternal 28,112,400; chr15 maternal 46,844,785; chr17 maternal 9,027,454; chr22 paternal
           35,690,655), each with 23 to 62 HiFi reads running through and none holding the 4 or 5 tracked k-mers.
           An alignment running through a place does not show that the read holds the assembly's sequence
           there; counting the k-mers in the reads does. At a wrong join the opposite problem appears: the
           reads' primary alignments go to the other haplotype, secondary records carry no sequence, and only
           alignment (no read through, several clipped) can speak.
           Along the way the test with simulated long reads found a crash on releases without any run of N, and
           the first full run stalled because it read a whole chromosome's track once per place.
Decision:  D1 settled: short reads detect, long reads adjudicate by k-mer content. With `--long-read-bam`, a run
           that long reads do not carry is a BASE_ERROR (reads run through) or a FALSE_JOIN (none does); a run
           they carry is no issue; a short copy-number loss on carried k-mers is dropped. Without long reads,
           whole-copy calls under 5 windows are low confidence.
Next:      the four unlisted places deserve a look by someone with the ONT and Element data.

### 2026-11-13 (j): rung 4, real assemblies; the classifier revised

Question:  on assemblies with real errors, are the copy-number issues real and how many discordances are found
           (GUIDE §8 rung 4; D4, D6)?
Run:       `scripts/run_hg002_assemblies.sh` on seven HG002 assemblies with the 36x read database, each scanned
           with the truth genome's k-mers as a second read set; `scripts/eval_truth.py`, `scripts/scoreboard.py`.
           Record: `experiments/2026-11-13_hg002_assemblies/`. Outputs: DATA/work/hg002_assemblies.
Seen:      With the classifier of entries h and i: high-confidence calls about 85% confirmed, losses found
           123 of 169 in `hprc_year2`, but collapses 8 of 64. The missed collapses were long stretches of
           duplicated sequence with a clean whole-copy pattern (4- and 6-copy k-mers at 4.9 and 7.3, 1- and 2-copy
           k-mers concordant), called ambiguous because each k-mer class was allowed the full regional wander of
           depth by itself. Three changes, each checked on the truth and planted genomes before the next:
           (1) one regional depth factor shared by all classes of a segment, each class within 4% of it:
           collapses found 17 of 64, truth genome 15 whole-copy segments; but four large planted events were
           lost, because (2) minor classes are mixtures and vetoed the fit: classes holding under 10% of the
           k-mers no longer veto; (3) the regional factor is bounded at two standard deviations, after a +1 call
           on the truth genome that needed a factor of 0.79.
           Then, in the scores by size: 18 high-confidence surplus-copy calls under 5 kb, 1 confirmed; the same
           family as 13 on the truth genome (entry k). Calls on fewer than 5 windows are now low confidence.
           And with the reads on `hprc_year2`: missing alleles 96% true when strong; one-sided missing sequence
           46% when strong, 9% when weak (planted genome: 51% and 10%). Weak one-sided clusters are no longer
           issues, which also took the truth genome's issue list from 1,297 rows to 404.
           Final, seven assemblies: 4,117 of 4,160 high-confidence copy-number issues confirmed (99%; 737 of 770
           without `verkko_2023`, which holds 150 Mb of surplus sequence). Surplus copies of 5 kb and more: 72%
           to 80% of assessable discordances found. Collapses: 26% to 50%. Under 5 kb: 10% to 20%.
           Depth anomalies: 19,446, 6% confirmed. Short-read k-mer QV sits 2 to 6 below the QV by truth k-mers.
           The numbers in entries h and i are those of the classifier before this entry; the final ones are in
           the experiment folders: truth genome 2,792 segments (2,630 depth anomalies, 20 whole-copy, 142
           ambiguous), 404 issue rows, 49 high-confidence at their locus, litmus 23 of 20,237 (18 high);
           planted genome 147 signals found (119 high), 13 of 656 issues unexplained.
Decision:  D4 and D6 as in `docs/method.md` §3. Thresholds were changed on rungs 2 to 4 together, so none of
           these rungs is a held-out test any more; rung 6 (a second sample, frozen thresholds) is the open one.
           (That test was run the next day, on NA19240: entries p and s.)
Next:      seven more assemblies are downloaded and listed; they wait for a disk that writes. (Six were run
           later the same evening, entry o; the numbers above are those of the first seven, and the scoreboard
           in the experiment folder holds all thirteen.)

### 2026-11-13 (i): rung 2, errors planted in the truth genome

Note:      the counts below are those of the classifier and reporting rules at that hour. Entry j changed both.
           The final counts are in `experiments/2026-11-13_planted_errors/` (147 signals found, 656 issues).
Question:  is each planted error found as the right class, and does anything else change (GUIDE §8 rung 2; D1, D4, D6)?
Run:       `scripts/plant_errors.py` (111 edits in hg002v1.1: 39 tandem duplications, 8 extra sequences, 36 deletions,
           12 swaps, 12 inversions, 4 inversions through N), `scripts/run_sample.py`, `scripts/eval_planted.py`.
           Record: `experiments/2026-11-13_planted_errors/`. Outputs: DATA/work/plant.
Seen:      Copy-number changes are found from 5 kb in unique sequence (duplication 4/4, deletion 4/4) and from 10 kb
           in duplications (4/4); below 5 kb mostly not (2 kb: 3/16). Identical extra sequences: 8/8 at 10 and 50 kb,
           partners 4/4. Missing haplotype after a deletion or a swap: 34/40 from 5 kb, 5/8 at 2 kb. Inversions:
           0/32 breakpoints as MISSING_JUNCTION (the homologue still holds the true adjacency); the 24 joined
           without N are all runs of unsupported assembly k-mers, the 8 joined through N have no k-mers to show.
           119 of 146 found signals are high-confidence. Of 1,582 issues, 276 belong
           to an edit, 1,269 are those of the unedited genome, 10 are at a homologous stretch, 27 are new (21
           MISSING_SEQUENCE, mostly 3 to 5 reads).
           Two defects showed up and were fixed here. (1) An HMM state that involves only k-mers of multiplicity 3 or
           more costs nothing in windows without such k-mers, so segments ran on (to the end of a simulated
           sequence) and were then labelled irregular; segments now use only windows that test their state. (2) An
           identical extra copy has no k-mers of its own, so "one copy too many" and "two thirds of the depth" could
           not be told apart (1/8 found at 10 kb).
Decision:  D4, D6. Added a step to S1: locate each candidate segment's k-mers in the release (`cpk locate
           --regions`), list the other copies, and test the k-mers unique to those copies at one copy against the
           segment's depth ratio (`copyproof/copies.py`). Ambiguous segments at least 8 windows long whose
           whole-copy reading changes depth 1.5-fold or more are reported at low confidence: the unedited genome has
           one such segment in 5.3 Gb, against 97 short dips to half depth. An unsupported run alone stays
           unreported (see h); inversions need the long-read audit.
Next:      rung 4. A planted genome with long reads simulated or re-aligned would test FALSE_JOIN; not done.
           (Done a day later with real reads aligned again: entry q.)

### 2026-11-13 (h): rung 3, the truth genome itself

Note:      the counts below are those of the classifier and reporting rules at that hour. Entries j and k
           changed both. The final counts are in `experiments/2026-11-13_hg002_q100_false_flags/` (2,792
           segments, 404 issue rows, 49 high-confidence at their locus; the litmus is unchanged at 0.11%).
Question:  what does copyproof flag on T2T-HG002 Q100 v1.1 with HG002's own reads (GUIDE §8 rung 3; D1, D5, D10)?
Run:       `scripts/run_hg002_truth.sh`, `experiments/2026-11-13_hg002_q100_false_flags/analyse.py`,
           `scripts/sd_boundary_litmus.py`. 711.6 million 151-bp reads (107.5 Gb). Record in that folder.
Seen:      Depth per copy 13.39. Assessable 88.7% (679 Mb masked: 94.8% of centromere and HSat annotation, 38.5% of
           segmental duplications). The HMM finds 2,791 departing segments; 2,597 are depth anomalies (one factor
           for all k-mer classes, or no regular pattern), 14 are whole-copy changes, 180 are ambiguous (27 of
           those reported). Litmus (denominator of entry g): 23 of 20,237 anchored duplication boundaries have an
           issue within 5 kb (0.11%; 19 high-confidence, 0.09%); target was 2%.
           Issues: 1,287, of which 94 high-confidence at their locus (47 per haplotype genome) and 646
           low-confidence, 558 of those MISSING_SEQUENCE with 3 to 5 reads. Read k-mers absent from the assembly:
           9.4 million solid ones, nearly all artefacts in four groups (adapter read-through; slippage at short
           tandem repeats; errors of k-mers read thousands of times, the mitochondrion alone giving 135,000;
           one-strand miscalls, mostly T to G, 13,102 loci). Assembly k-mers with no read: 13,044 runs, 1,123 with
           covered flanks; of 160 sampled, HiFi reads run through 157 and the other 3 sit at rDNA gap edges.
           Short-read k-mer QV of this genome: 52.4, a floor of the chemistry, not the assembly's quality.
           High-confidence copy-number calls include the first 60 kb of chr1 on both haplotypes (3 to 4 copies more
           in the reads), eight 1 to 2 kb stretches whose haplotype-specific k-mers have no reads (two on known
           Q100 issues), and immunoglobulin loci (V(D)J deletion in the B-cell line: source discordance).
Decision:  D5 result 0.11%. D10: no change of course. D1: short reads detect; long reads are needed to say whether an
           unsupported run is a false join, and help for single-class half-depth segments. D6: a depth factor is
           not an issue; ambiguous segments need support; unsupported runs need a second witness.
Next:      rung 2.

### 2026-11-13 (g): D5 fixed before the 35x truth run

Question:  what is the denominator of the false-flag litmus (GUIDE §8 rung 3, decision D5)?
Run:       none; written down at 14:20, before the 35x database existed.
Seen:      n/a.
Decision:  D5. Denominator: every boundary of an interval of `hg002v1.1.SDs.013025.merged.bed` whose flank outside the
           duplication is at least 5 kb long (to the next duplication or the sequence end) and free of centromere and
           HSat annotation within 5 kb. A boundary is flagged when a copyproof issue of class COLLAPSE, EXPANSION,
           UNSUPPORTED, MISSING_JUNCTION, MISSING_HAPLOTYPE or FALSE_JOIN lies within 5 kb. Report the rate for all
           issues and for high-confidence issues, next to issues per haplotype genome outside the mask.
           Script: `scripts/sd_boundary_litmus.py`.
Next:      run rung 3 at 35x.


### 2026-11-13 (f): S2, read k-mers the release lacks

Question:  can solid read k-mers absent from the release be localised and classified (GUIDE §3 S2)?
Run:       `copyproof/missing.py` on the 18x and 36x truth runs; `cpk recruit`, `cpk locate`.
Seen:      At 18x the first version recruited 2 million reads (its cap) for a genome that should lack nothing.
           Filters were added one at a time, each after looking at what the reads were: adapter trimming at
           counting (TruSeq read-through makes solid k-mers found in no genome); a low-complexity filter (tandem
           units of 1 to 6 bases); a ceiling of 3.5 copies; a set of high-copy k-mers by multiplicity and by read
           depth, applied per read and per hit; a strand-balance filter. Clustering by shared k-mers made one giant
           cluster; clustering by anchored locus does not. Recruiting took 708 s until untracked k-mers were
           skipped; classifying 622,000 reads took over 10 minutes until chaining was vectorised (33 s).
           A junction is placed within 12 bp on simulated data.
Decision:  S2 as in `docs/method.md` §4. MISSING_HAPLOTYPE (three or more missing alleles within 10 kb) is the most
           useful class: it finds a missing or overwritten haplotype with no parents.
Next:      rung 3 at 36x.

### 2026-11-13 (e): S1, how to tell a copy change from a depth change

Question:  what model of read counts along the assembly separates copy-number errors from uneven depth (D4)?
Run:       truth genome at 18x (one read file), then 36x; `copyproof/cn.py`.
Seen:      A first HMM with whole-copy states only called about 1,000, then 2,400 segments on the truth genome:
           every dip in depth became a lost copy. Changes, in order: a mask (too few usable k-mers, more than half
           above the multiplicity cap, GC factor below 0.5, short islands); multiplicative states for depth
           factors next to the additive ones; a classifier on the summed evidence of each segment with a test of
           adequacy, not only a likelihood ratio; the commonest class must itself follow a whole-copy change;
           "absent" only for 20 windows or more (263 short dropouts otherwise, 1 true); which classes follow a
           change is part of the hidden state (drawn per window, it taxed every window and hid collapses among
           high-copy k-mers). A local depth baseline from low-copy k-mers was tried and switched off: the copies of
           a multi-copy k-mer lie elsewhere, so local depth does not predict them (truth-genome calls 17 to 40).
Decision:  D4: an HMM over 1-kb windows with a concordant state, additive states (+1, -1, +2, each for k-mers of
           multiplicity at least 1, 2, 3 or 5) and depth factors (0, 0.25, 0.5, 0.75, 1.25, 1.5, 2), Student-t
           noise, then classification per segment into CN, CN_OR_DEPTH, DEPTH, ABSENT, IRREGULAR. Emissions are not
           negative binomial as the guide leaned: the variance is derived (entry d).
Next:      S2.

### 2026-11-13 (d): the count model

Question:  what is the variance of a window's summed k-mer counts (GUIDE §3 S1)?
Run:       truth genome, 18x and 36x; `copyproof/model.py`, `copyproof/windows.py`.
Seen:      Fitting one overdispersion per class gave 118: neighbouring k-mers share reads, so their counts are far
           from independent, and heterozygous k-mers come in runs of up to 31. The covariance of two k-mers d bases
           apart is depth x (1 - d / R), R being the k-mer starts per read; summing it over all pairs of a class in
           a window gives a weight V, and Var(S) = omega x mu x V + (psi x E)^2 fits with omega between 1.0 and 4.4
           and psi = 0.077 at 36x. Checked against brute force. Copies per class come out within 2.5% of the
           multiplicity for classes 1 to 8.
           Depth per copy needs no peak fitting: it is half the median count of k-mers found once in each
           haplotype (13.39 at 36x). X and Y in a male need no special case.
Decision:  the model of `docs/method.md` §2.
Next:      segmentation.

### 2026-11-13 (c): rung 0, a synthetic trio

Question:  does each channel find a planted error of its class in simulated data (GUIDE §8 rung 0)?
Run:       `copyproof/sim.py`, `tests/` (30 tests, about a minute).
Seen:      A 600-kb diploid genome with a collapsed duplication, a false duplication, a dropped stretch, a
           haplotype switch and 10 kb of foreign sequence: each is found once, on the right haplotype, as the right
           class; the homologous stretch is a `partner`; the de novo junction is placed within 12 bp and is absent
           from both parents; nothing else is flagged. A second case holds an identical extra sequence.
Decision:  the tests are the regression suite; every later change ran against them.
Next:      real data.

### 2026-11-13 (b): what it costs

Question:  wall time and memory on whole genomes (D9).
Run:       Apple M1 Max, 10 cores, 64 GB. `experiments/2026-11-13_engine_validation/`.
Seen:      Counting 107.5 Gb of gzipped FASTQ: 304 s, 20.4 GB. Counting a 30x CRAM through samtools: 29 minutes
           (decoding, with the machine busy). Scanning a 6-Gb diploid release: 93 s, 17.3 GB. Windows: 71 s with 6
           processes. Copy number with the copy lookup: about 3.5 minutes. Missing k-mers: one more pass over the
           reads. A k-mer database is 6.5 to 7 GB and serves every assembly of the sample.
Decision:  D9: counting dominates, as expected, and is cheap enough that region-mode scans are not needed.
Next:      none.

### 2026-11-13 (a): one k-mer engine instead of existing counters

Question:  reuse meryl, KMC or jellyfish, or build (D8)?
Run:       `engine/` (Rust), `tests/test_engine.py`, `experiments/2026-11-13_engine_validation/`.
Seen:      Every channel needs the same two products: read counts of the assembly's k-mers laid out along the
           assembly, and the read k-mers the assembly lacks. A counter of all k-mers needs hundreds of GB of
           temporary space per read set. Tracking the k-mers whose hash is a multiple of 4 and counting those
           exactly gives the same information (neighbouring k-mers share reads) in one pass and 20 GB. Against
           jellyfish 2.3.1 on 4 million reads: 13,534,076 of 13,534,369 tracked k-mers identical, 292 higher by one
           (Bloom filter), one saturated; none missing, none lower.
Decision:  D8: build (`cpk count | scan | recruit | locate`). D2: k = 31 throughout; other sizes were not tried.
Next:      the count model.

### 2026-11-13: guide written, resources checked

Question: what public data can drive each rung of the validation ladder (GUIDE §8)?

Run: listings of the HPRC data tables, the HGSVC3 release and working folders, and the HG002 Q100 folders. Every entry in resources.tsv was checked on this date.

Seen:

1. HPRC publishes, for every R2 haplotype, gap BEDs, segdup and censat annotation, chains to GRCh38 and CHM13, Flagger's read-to-assembly BAMs and high-clip bigWigs. S3–S5 can therefore start on HPRC R2 without realignment.
2. A known-issue BED exists for hg002v1.1 (3,172 intervals, 38.8 Mb). Rung 3 can separate known truth issues from new flags.
3. hg002v1.2 is public, with a chain to v1.1. Loci fixed between the versions are natural positives.
4. 24 HG002 assemblies have GQC outputs against v1.1, giving real, labelled errors for rung 4.
5. HGSVC3 ships unassigned and disconnected sequence in separate files. S1 needs to count it, or it will report false collapses.

Decision: none yet; these are inputs for M0–M2.

Next: M0. Re-run both case drivers and checksum the HG002 inputs.
