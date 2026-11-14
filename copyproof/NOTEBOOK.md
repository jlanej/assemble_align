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
| D1 | Read source for copy number | short reads detect; long reads adjudicate by k-mer content (`--long-read-bam`) | h, k |
| D2 | k-mer size | 31 throughout; other sizes not tried | a |
| D3 | Coordinate frame | the release's own sequences; nothing is projected to a reference | a |
| D4 | Segmentation | HMM over 1-kb windows with additive and depth-factor states, then classification per segment with one shared regional factor | e, i, j |
| D5 | Negative-control denominator | fixed in g before the run; result 0.11% | g, h |
| D6 | Thresholds and score | rules with log-likelihood margins; changed on rungs 1 to 4, so not frozen: rung 6 is open. Every value and its origin: `docs/thresholds.md` | e, h, i, j, l, m, n |
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
Entries m, n and o were written within the hour of their experiments.
The experiment folders are regenerated from the run directories by the scripts named in each, and are the
reference where a number here and a number there differ.

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
