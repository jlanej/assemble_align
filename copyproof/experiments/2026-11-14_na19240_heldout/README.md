# Rung 6: NA19240, the held-out test

**Question.** With every threshold fixed beforehand, do the results hold on a second person
(GUIDE §8 rung 6)? Rungs 1 to 4 adjusted thresholds on the data they report. This test was
written down first, committed and pushed, and run afterwards.

In the guide's ladder the comparison itself is rung 5, agreement without truth, for one of the
four people who have both an HPRC and an HGSVC3 assembly. It was run under the condition that the
guide sets before rung 6, the scan of whole releases: every threshold frozen. The notebook and the
other documents call it rung 6, and so does this record. The scan itself was not run.

**Answer.** All seven comparisons that were declared hold. That is a weaker statement than it
sounds, and the rest of this record says why: there is no truth genome for this person, the
stand-in for one moved by 19 points between two pairs of assemblies on HG002, and it cannot judge
collapses.

| Score | Predicted from HG002 | NA19240 | p, one-sided | Holds |
|---|---|---|---|---|
| 1. High-confidence copy-number issues confirmed by both other assemblies | 74 of 119 | 46 of 81 | 0.27 | yes |
| 2. Discordances of 5 kb and more found, the others hold more (collapse) | 1 of 3 | 2 of 14 | 0.46 | yes |
| 2. The same, the others hold fewer (surplus copy) | 117 of 127 | 64 of 76 | 0.066 | yes |
| 3. Planted false duplications found, from 5 kb | 46 of 47 | 47 of 48 | 0.76 | yes |
| 3. Planted deletions found as a collapse on the homologue, from 5 kb | 25 of 28 | 28 of 31 | 0.71 | yes |
| 3. Planted deleted or overwritten stretches found as a missing haplotype, from 5 kb | 34 of 40 | 41 of 43 | 0.98 | yes |
| 3. Planted swaps found as a phase switch | 7 of 12 | 8 of 12 | 0.80 | yes |

(`verdicts.tsv`. A score holds unless NA19240's share is below the predicted one by a one-sided
Fisher exact test at 5%; seven comparisons, no correction for their number.)

## What was fixed, and when

Everything is in `NOTEBOOK.md`, entry p, written at 07:34 on 2026-11-14 (CST), before any read of
NA19240 was counted and before any program had opened her assemblies. Later additions to the
entry carry their times. `prereg.sha256` holds the checksum of each wording and of every script,
with the time each line was added.

| Time | What |
|---|---|
| 07:34 | Entry p: sample, release, evidence, the three scores, the rule for the verdict |
| 07:39 | One assertion added to the test of the scoring script: the calls do not change when label sets ride along |
| 07:50 | The driver waits for a download to be complete before it counts |
| 07:53 | The predictions, after the calibration on HG002 (`predictions.tsv`) |
| 08:00 | Two points from a review: the hapmer audit is not an independent test of a trio-phased release; counts are reported, not only shares |
| 08:05 | One run added on request: the hapmer audit on the Strand-seq-phased Verkko assembly |
| 08:25 | Four lines of the entry reworded; no rule, number or time changed |
| 10:22 | `tables.py`: its half that reads NA19240's runs was run for the first time, on HG002's run under NA19240's file names, and two faults in reading files were fixed |
| 10:27 | Committed on request as `a82fdbe`, merged with the branch as `8f5b972` |
| 10:29 | Pushed. GitHub's activity record gives 16:29:27 UTC |
| 10:30 | The scored run starts (10:30:45); it had been held for the push |

The code under test is commit `1465496` with every default it carries. The drivers refuse to run
when `src/`, `engine/src` or `scripts/run_sample.py` differ from it; they never did. Of the two
values that were not the best on the old data, the shipped ones were tested (`min_seg_class` 20,
`margin` 3). No threshold was changed after the first result, and no code under `src/` or
`engine/`.

## Inputs

- **The release under test:** `NA19240_pat_hprc_r2_v1.0.1.fa.gz` and
  `NA19240_mat_hprc_r2_v1.0.1.fa.gz` of HPRC release 2 (hifiasm 0.19.9, trio): 64 and 65
  sequences, 6,063,484,360 bases. copyproof had not been run on this person.
- **Reads:** the NYGC 30x CRAMs of NA19240 (727,171,320 reads), her mother NA19238 (939,885,340)
  and her father NA19239 (723,690,456), all of 150 bases, each counted once.
- **Two other assemblies of the same person**, from HGSVC3, made at another centre from other
  reads and phased with Strand-seq, without the parents: a Verkko assembly (1,120 sequences in
  seven files, 6,106,495,734 bases) and a hifiasm assembly (248 sequences in two files,
  7,171,149,766 bases).
- **HPRC's own material for this release:** the annotation of segmental duplications and
  satellites, the Flagger and NucFlag calls, and the HiFi-to-assembly alignment, read remotely.

Sizes and counts are in `kmer_databases.tsv` and `summary_short.json`. The manifest is
`scripts/manifests/05_na19240_heldout.tsv`. Every assembly and CRAM was checked against its
published MD5 (`inputs.md5`, `check_inputs.sh`). That check was needed: one file of the hifiasm
assembly arrived from the EBI server with the right size and 198 wrong bytes, an error message of
the server in place of data. Five of six later requests for that block returned the right bytes;
the block was replaced and the file then matched its checksum.

## The stand-in for a truth genome, and what HG002 says about it

No truth genome exists for NA19240. The exact k-mer counts of the two HGSVC3 assemblies label
every tracked k-mer of the release as a truth genome's do in rung 4. A k-mer is contradicted *by
both* when the two assemblies hold the same number of copies of it and that number differs from
the release's. That rule, `both`, is the declared one. Each count is also given for the Verkko
assembly alone, the hifiasm assembly alone, and either.

This is a proxy, not truth. It was applied to HG002 first, where the truth is known: the HPRC
year-2 assembly, judged by a Verkko assembly (`v2_trio`) and a hifiasm assembly from other reads
(`hprc_year1`). Of its 119 high-confidence copy-number issues (`hg002_proxy_against_truth.tsv`,
`score1_counts.tsv`):

| Judged by | Confirmed | Of those, not confirmed by the truth |
|---|---|---|
| the truth genome | 114 | |
| both, with the same number of copies | 74 | 0 |
| the Verkko assembly alone | 101 | 2 |
| the hifiasm assembly alone | 82 | 0 |
| either | 103 | 2 |

So `both` confirmed nothing that the truth did not, and missed 40 that the truth confirmed. A
single assembly came nearer to the truth's count and confirmed 2 of the 5 issues that the truth
does not. Collapses are where the proxy fails: of 24 high-confidence collapse issues the truth
confirmed 22, the Verkko assembly 9, the hifiasm assembly 2 and `both` none.

A second pair of HG002 assemblies (`verkko_asm_thic`, `hprc_year2_v2_thic`) gave `both` 97 of 119.
The proxy thus moved by 19 points with the assemblies chosen, before NA19240 was run. The
declared prediction was the lower of the two.

## Run

```bash
D=/path/to/data
scripts/fetch_manifest.sh scripts/manifests/05_na19240_heldout.tsv $D
experiments/2026-11-14_na19240_heldout/check_inputs.sh $D
experiments/2026-11-14_na19240_heldout/calibrate_hg002.sh $D                 # HG002, then the predictions
scripts/run_na19240_heldout.sh $D 10 count refs short score_short            # the scored run
scripts/run_na19240_heldout.sh $D 10 long score_long                         # with HPRC's HiFi alignment
experiments/2026-11-14_na19240_heldout/plant.sh $D 10                        # score 3
scripts/run_na19240_heldout.sh $D 10 verkko_labels verkko score_verkko       # the added run
.venv/bin/python scripts/dev/cache_segments.py $D/work/rung6/NA19240_short
.venv/bin/python experiments/2026-11-14_na19240_heldout/alt_settings.py $D/work/rung6/NA19240_short \
    $D/work/rung6/NA19240_short.eval experiments/2026-11-14_na19240_heldout/alt_settings_na19240.tsv
.venv/bin/python experiments/2026-11-14_na19240_heldout/tables.py $D experiments/2026-11-14_na19240_heldout
.venv/bin/python experiments/2026-11-14_na19240_heldout/after.py $D experiments/2026-11-14_na19240_heldout
```

The scored run (`short`) is the child's reads and both parents, without long reads, because
neither the calibration nor rung 4 has long reads. It took 1,167 s, of which 826 s were the pass
over the reads for missing sequence (`runs.tsv`). The run with long reads took 1,251 s on top of
the shared scan, 1,137 s of that for fetching and reading alignments; the planted genome 1,011 s;
the Verkko assembly 329 s (`summary_*.json`). The three read counts took 2,721 s, 2,450 s and
758 s for 109, 141 and 109 Gb of reads (`kmer_databases.tsv`). The first two ran while a CRAM was
being downloaded to the same drive, and other sessions used the machine all morning; the third
ran with no download. Which load cost the time was not measured, so these are not costs of the
tool.

## Score 1: are the high-confidence copy-number issues real?

The scored run reported 81 high-confidence copy-number issues: 68 surplus copies (`EXPANSION`),
12 collapses and 1 `UNSUPPORTED` (`issues_by_class.tsv`). By the other assemblies
(`score1_counts.tsv`):

| Judged by | Collapses (12) | Surplus copies and unsupported sequence (69) | All (81) | Not confirmed |
|---|---|---|---|---|
| both, with the same number of copies | 2 | 44 | 46 | 35 |
| the Verkko assembly alone | 11 | 64 | 75 | 6 |
| the hifiasm assembly alone | 9 | 44 | 53 | 28 |
| either | 11 | 64 | 75 | 6 |

The 35 that `both` does not confirm are listed one by one in `unconfirmed_high.tsv` and sorted in
`after_confirmed_by.tsv`:

- **22 are confirmed by the Verkko assembly only** (20 surplus copies, 2 collapses). The hifiasm
  assembly's discordances cover none of them, not even in part.
- **7 are confirmed by each assembly, with different numbers of copies.** All 7 are collapses.
- **6 are confirmed by neither** (5 surplus copies, 1 collapse). Two of them, of 73 kb and 442 kb,
  are covered for 41% and 35% by the Verkko assembly's discordances and lie at Flagger and at
  NucFlag calls. The other four are of 5 to 92 kb, rest on 5 to 10 windows, and nothing supports
  them.
- None is confirmed by the hifiasm assembly only.

**The declared comparison holds, and says little.** 46 of 81 (57%) against the predicted 74 of
119 (62%) is not lower at 5% (p = 0.27). Against the second HG002 pair, 97 of 119 (82%), it is
lower (p = 0.0001; `after_against_hg002.tsv`). The difference between the two readings is the
second judge. On NA19240 the hifiasm assembly alone confirms 44 of the 69 where the Verkko
assembly confirms 64. On HG002, `hprc_year1`, a hifiasm assembly, confirmed 80 of 95, and
`hprc_year2_v2_thic` 92 of 95. The second pair seems to hold no hifiasm judge at all:
`verkko_asm_thic` is a Verkko assembly, and the sequences of `hprc_year2_v2_thic` carry names of
Verkko's kind (`mat-0000001`), though rung 4's record does not name its assembler.

**The hifiasm assembly of HGSVC3 is a poor judge by itself.** It holds 7.17 Gb, about 1.1 Gb more
than the Verkko assembly, in about as many distinct k-mers (`kmer_databases.tsv`). Alone it
reports 2,818 assessable stretches of 5 kb and more where it holds more copies than the release,
against the Verkko assembly's 129 (`na19240_short_cn_scores.tsv`). No issue covers 2,813 of them,
and at 2,777 of those, 1.38 Gb together, the read depth is nearer to what the release holds than
to what the hifiasm assembly holds (`after_reads_as_referee.tsv`). The surplus is the hifiasm
assembly's. The rule `both` was not misled by it, since the Verkko assembly must hold the same
count. But an assembly that carries surplus copies of its own cannot contradict a surplus copy in
the release, and that is what lowers `both` here.

**By the Verkko assembly alone, 75 of 81 (93%) are confirmed.** Its counterparts on HG002
confirmed 101 and 98 of 119 (85% and 82%), where the truth confirmed 114 (96%). So neither
measure shows a drop from HG002 to NA19240. Neither is the true share: on HG002 one assembly
confirmed 2 of the 5 issues that the truth does not, and failed to confirm 15 of the 114 that it
does.

## Score 2: are the discordances found?

By `both`, the release differs from the two other assemblies in 90 assessable stretches of 5 kb
and more: 76 where it holds more copies and 14 where it holds fewer
(`na19240_short_cn_scores.tsv`).

**Surplus copies: 64 of 76 found, against 117 of 127 predicted (p = 0.066).** This holds, but
narrowly; a real drop of some points cannot be excluded. `after_missed_summary.tsv` and
`after_missed_discordances.tsv` say what the copy-number channel holds at the 12 that no issue
covers:

| What the channel holds there | Stretches | Read depth over the release's count | The others' count over the release's |
|---|---|---|---|
| a segment judged a depth change, not a copy change | 3 | 0.65 | 0.67 |
| an ambiguous segment, which is not reported | 4 | 0.86 | 0.86 |
| an irregular segment | 3 | 0.17 | 0.64 |
| no segment | 2 | 1.00 | 0.82 |

(Medians.) In 7 of the 12 the reads show the drop that the other assemblies predict, and the
classifier declined to call it. In 3 the reads hold almost none of the k-mers. In 2 the reads
agree with the release. At the 57 found by a whole-copy segment the two ratios are 0.674 and
0.677.

**Collapses: 2 of 14 found, against 1 of 3 predicted.** This comparison cannot show anything, as
entry p said when three discordances were all that HG002 offered. The 12 not found do not look
like missed collapses. At 11 of them the channel has no segment. At 10 of those the reads stand at
1.04 to 1.15 times the release's count, where the other assemblies stand at 1.4 to 2.0 times; at
the eleventh the reads stand at 0.76. Ten of the twelve are five pairs: the same stretch on both
haplotypes by the alignment of the two, discordant on both and found on neither, among them
390 kb in two pieces. The reads that copyproof was given agree with the release there. Whether
both HGSVC3 assemblies carry the same surplus, or their DNA differs from the DNA that NYGC
sequenced, was not examined.

**By the Verkko assembly alone, 102 of 146 surplus-copy stretches are found (70%)**, against 129
of 153 (84%) and 129 of 163 (79%) on HG002, and that is lower (p = 0.002 and 0.04;
`after_against_hg002.tsv`). It was not a declared comparison, and a single judge's discordances
include its own shortfalls: on HG002, 57 of the 228 of `v2_trio` were not discordances by the
truth (`hg002_proxy_against_truth.tsv`). The reads can referee: of the 44 not found, the read
depth is nearer to the release at 23 and nearer to the Verkko assembly at 21
(`after_reads_as_referee.tsv`). Counting only stretches where the reads side with the Verkko
assembly, 95 of 116 are found. That sorting was made after the result was seen and uses the read
depth, which the calls use too, so it is a description and not a score.

## Score 3: planted errors

115 of 116 edits of rung 2's design were placed in the release (no place was found for one 200-kb
deletion in a segmental duplication), on sequences of at least 12 Mb, outside satellites, outside
HPRC's Flagger and NucFlag calls and outside the mask of the scored run (`planted_edits.tsv`).
The planted genome was run like the unedited one and scored by `scripts/eval_planted.py`,
unchanged (`planted_sensitivity.tsv`). All four declared counts hold and none is below HG002's
share.

What the short reads do not find is what they did not find on HG002: no deletion as a false join
(0 of 39), no inversion breakpoint as a missing junction (0 of 24), and no deletion of 2 kb as a
collapse (0 of 8). The 63 planted false adjacencies of deletions and inversions all show as runs
of assembly k-mers that no read holds; the 8 joins through N do not (`planted_breakpoints.tsv`).
Of 2,700 issues on the planted genome, 249 lie at edits and 2,434 are also on the unedited
release. 17 are new, 4 of them near the stretch homologous to an edit; on HG002 it was 13 and 1
(`after_planted_totals.tsv`; `planted_new_elsewhere.tsv` lists the other 13).

## Reported without a prediction

**Issues by class** (`issues_by_class.tsv`, beside HG01943 of rung 1). The scored run reports
2,438 issues, the run with long reads 2,970.

**Missing alleles.** Of 1,773 strong missing-allele calls, the Verkko assembly holds the k-mers of
1,259, the hifiasm assembly of 1,211, both of 1,075 and either of 1,395 (79%;
`na19240_short_missing_clusters.tsv`). On HG002 either held 913 of 996 (92%) and the truth 957.
The difference lies in one kind of call (`after_missing_alleles_by_parent.tsv`; finer in
`after_missing_alleles.tsv`):

| | Calls | Held by one or both other assemblies | Median read count |
|---|---|---|---|
| NA19240, k-mers in a parent's reads | 851 | 774 (91%) | 14 |
| NA19240, k-mers in neither parent's reads | 922 | 621 (67%) | 8 |
| HG002, k-mers in a parent's reads | 925 | 857 (93%) | 13 |
| HG002, k-mers in neither parent's reads | 71 | 56 (79%) | 8.5 |

The depth of one copy is 13.6 on NA19240 and 13.4 on HG002. For alleles that a parent's reads
carry, the two samples agree. NA19240 has thirteen times as many calls in neither parent's reads,
at about six tenths of the depth of one copy, and the other assemblies hold two thirds of them.
What they are was not examined; variation within the cell line and an effect of the child's
library are both open.

**Assembly k-mers that no read holds.** 203,153 of the release's tracked k-mer positions are in
no read. The Verkko assembly lacks 67,926 of them, the hifiasm assembly 67,631, both 65,309
(`na19240_short_unread_kmers.tsv`). On HG002 both lacked 105,427 of 281,254 and the truth 141,920.

**Hapmers on the release.** 88 sequences assessed and 41 not; 88 pieces between gaps, 85
consistent with their sequence, 3 not assessable and none on the other parent; 36 switch
segments over 559 kb, 18 of them high-confidence (`runs.tsv`). The release was phased with these
parents, so few switches do not show that its phasing is right.

**HPRC's Flagger and NucFlag calls.** `issues_at_published_calls.tsv` and
`published_calls_at_issues.tsv` are of the run with long reads; the `after_*` tables are of the
scored run.

- The issues lie at NucFlag calls far more often than issues placed at random on the same
  sequences would (`after_issues_at_calls.tsv`, 200 placements): all 18 high-confidence phase
  switches (0.9 expected, at most 3), 52 of 76 high-confidence missing haplotypes (1.2), 59 of 76
  high-confidence base errors (1.0), 33 of 68 high-confidence surplus copies (6.9), 4 of 12
  high-confidence collapses (0.1). Nearly all of these lie at calls labelled `MISJOIN`
  (`after_issues_by_call_label.tsv`).
- Flagger's calls touch 18 of the 68 high-confidence surplus copies (4.1 expected) and 1 of the 12
  collapses (1.0 expected).
- The other way round, 80 of Flagger's 90 `Col` and `Dup` calls lie mostly in sequence that
  copyproof does not assess; of the 10 that do not, 5 are at a copy-number issue. Of its 59 `Err`
  calls, 37 are assessable and 29 of those are at a copy-number issue
  (`after_published_calls.tsv`).

**The run with long reads.** 25 joins through N have no linking read and 30 cannot be assessed.
Of 3,000 runs of unread assembly k-mers the long reads carry 505, do not carry 519, are unclear
at 1,973, and 3 could not be fetched (`runs.tsv`). The scores of this run as
`scripts/eval_assemblies.py` gives them (`na19240_long_cn_scores.tsv`: 76 high-confidence
copy-number issues, 42 confirmed by `both`, 70 by the Verkko assembly) leave issues out: the
script scores issues whose channel is exactly `S1`, and the long-read audit marks the 38 it looked
at `S1+S4`. `after_long_reads_summary.tsv` follows every copy-number issue of the scored run into
this one:

- With every channel counted there are 79 high-confidence copy-number issues; `both` confirms 44
  and the Verkko assembly 72 (`after_long_reads_score1.tsv`).
- Of the 69 high-confidence surplus copies, the audit looked at 5. At 3 the long reads carry the
  k-mers, and the issue was made low-confidence; the Verkko assembly confirms all 3, `both` 2. At
  2 they do not, and the issue stayed; both are confirmed.
- Of the 82 low-confidence ones, 32 were marked as carried by the long reads (the Verkko
  assembly confirms 1), 12 are no longer reported (it confirms none), and 1 was raised to
  high-confidence (confirmed by neither).

## The two candidate values

The segments of the scored run, classified again as the sensitivity record did
(`alt_settings_na19240.tsv`; without the lookup of other copies, so "as shipped" stands at 78
high-confidence issues and not 81):

| Setting | High-confidence | Confirmed by both | By the Verkko assembly | Surplus copies of 5 kb found, of 76 |
|---|---|---|---|---|
| as shipped (20, 3) | 78 | 45 | 72 | 64 |
| `min_seg_class` 10 | 77 | 44 | 72 | 66 |
| `margin` 2 | 84 | 47 | 78 | 65 |
| both | 84 | 49 | 79 | 67 |

With both values changed there are 6 more high-confidence issues, 7 more confirmed by the Verkko
assembly and 3 more surplus copies found. On HG002 the same change gave 22 more issues and 20
more confirmed by the truth (`alt_settings_hg002.tsv`). The direction is the same on the new
sample. This was declared as a second result; it changes neither the scored run nor the shipped
values.

## The added run: the hapmer audit where it is independent

The Verkko assembly was phased with Strand-seq, so an audit with the parents' reads tests phasing
that the parents did not help to make. Every file of it was run at the frozen thresholds with the
same reads; HPRC's two trio-phased haplotypes rode along as label sets, and each call was held
against them by the rule of entry p (`verkko_hapmers_summary.tsv`).

- **Sequences.** 177 of 1,120 sequences were assessed, 6.07 Gb. For all 103 assessed sequences of
  the two haplotype files, hapmers and HPRC's k-mers name the same parent. Of the 74 assessed
  among the other files, 24 agree, 20 disagree (7.6 Mb) and 30 are unclear
  (`verkko_hapmers_sequences.tsv`).
- **Switches.** 65 segments over 27.2 Mb. Of the 13 with 10 sites or more, HPRC's haplotypes
  confirm 7, contradict none and are unclear at 6. Of the 52 with fewer sites they confirm 16,
  contradict 8 and are unclear at 28.
- **The largest.** Three confirmed segments hold 24.8 of the 27.2 Mb
  (`verkko_hapmers_calls.tsv`): the last 12.4 Mb of a maternal sequence of the second haplotype
  file follow the father (8,164 sites against 29), and 12.4 Mb at the end of an unassigned
  sequence that is otherwise paternal follow the mother (8,769 against 17).
- **The other way round.** Of 121,547 bins of 50 kb, 660 follow the other parent by HPRC's
  k-mers; the audit flags 482 of them. Of the 542 bins it flags, 482 follow the other parent
  (`verkko_hapmers_other_parent_bins.tsv`).
- No piece between gaps was called misplaced.

HPRC's haplotypes were phased with the same parents' k-mers from other reads: another assembly,
not another pedigree.

## What was learned

1. **Frozen, the thresholds gave on a second person what they gave on HG002**, by every
   comparison that was declared. The planted errors, which need no proxy, are found at the same
   rates or higher.
2. **The share of high-confidence copy-number issues that are real is not pinned down for
   NA19240.** The two measures give 57% (both assemblies, same count) and 93% (the Verkko
   assembly alone); on HG002 they gave 62% and 85% where the truth gave 96%. Six of 81 issues are
   confirmed by neither assembly.
3. **Surplus copies of 5 kb and more are found at 84% by the declared measure, against 92% on
   HG002.** The difference is within what the test allows (p = 0.066). Seven of the twelve misses
   are segments that the classifier saw and did not call; the two candidate values recover three.
4. **Collapses cannot be assessed with this proxy.** Where the two other assemblies hold more
   than the release and no collapse is called, the reads agree with the release at 10 of 12.
5. **A second assembly is a judge only as far as it is free of the fault in question.** The
   HGSVC3 hifiasm assembly carries 1.38 Gb of surplus copies in its assessable stretches of 5 kb
   and more, and so confirms fewer surplus copies in the release. One assembly from another
   assembler was the more useful judge on both samples.
6. **Half of NA19240's strong missing-allele calls are in neither parent's reads**, at about six
   tenths of one copy's depth; on HG002 7% are. With parents at hand such calls can be told
   apart. The tool does not do so.
7. **On an assembly made from the same long reads, the long-read audit withdrew confidence from
   surplus copies that the Verkko assembly confirms** (3 of the 5 it looked at). The counts are
   small.
8. **The hapmer audit finds phase switches in an assembly phased without the parents**: 7 of its
   13 high-confidence segments are confirmed by trio-phased haplotypes and none is contradicted.
   Its low-confidence segments are contradicted at 8 of the 24 that get a verdict.

## What this test does not show

- That the calls are true. No truth genome exists for this person, and each statement above is
  agreement between assemblies or between assemblies and reads.
- Anything about collapses, beyond the planted ones.
- Anything beyond one more person. The release is of the same assembler and pipeline as HG01943
  and as HG002's `hprc_year2`, and the reads are of the same centre, read length and depth as
  HG01943's.
- Independence of the two judges from each other: both HGSVC3 assemblies were made from the same
  reads.
- That the scoring rule was neutral. It was written the same morning, before NA19240 was opened
  but after HG002 had been seen, and its prediction was the lower of two anchors.

## What went differently from the plan

- The pre-registration was to rest on local times only. It was committed and pushed on request
  before the scored run (see the table above); entry p records both statements.
- `tables.py` was changed at 10:22, after its checksum had been logged at 08:05 and before any
  run of NA19240: its untested half was tried on a stand-in and two faults were fixed. Its
  output for HG002, the predictions among it, stayed byte for byte the same.
- The scores of the run with long reads leave out the issues that the long-read audit marked
  (above). The scoring script was not changed; the complete counts are in the `after_*` tables.
- `after.py` and every `after_*` table were written after the scored run had been seen. They
  describe and decide nothing.
- HPRC's Flagger and NucFlag files were to stay unopened until the runs were done. Their first
  three lines were read to learn the format, and nothing else before the scored run.
- One comment line in each of `scripts/run_na19240_heldout.sh` and `hapmers_against_hprc.py` was
  corrected after the last run; `git diff a82fdbe` shows them.
- A second genome with a near-truth assembly, I002C, was considered and not used: its public
  short reads are MGI and its parents' short reads are not public (`resources.tsv`).

## Files

| File | What | Written by |
|---|---|---|
| `FROZEN_COMMIT`, `prereg.sha256`, `predictions.tsv` | the pre-registration | by hand; `tables.py` |
| `inputs.md5`, `check_inputs.sh` | published checksums of the inputs, and the check | by hand |
| `calibrate_hg002.sh` | the proxy applied to HG002 | |
| `hg002_cn_scores.tsv`, `hg002_proxy_against_truth.tsv`, `hg002_missing_clusters.tsv`, `hg002_unread_kmers.tsv`, `hg002_pair2_*.tsv`, `hg002_planted_with_parents.sensitivity.tsv` | the calibration | `scripts/eval_assemblies.py`, `scripts/eval_planted.py`; copied by `tables.py` |
| `verdicts.tsv`, `score1_counts.tsv` | the seven comparisons; score 1 as counts for every judge | `tables.py` |
| `na19240_short_*.tsv`, `na19240_long_*.tsv` | scores, missing-sequence clusters and unread k-mers of the two runs | `scripts/eval_assemblies.py`; copied by `tables.py` |
| `unconfirmed_high.tsv` | the 35 high-confidence issues that `both` does not confirm | `tables.py` |
| `issues_by_class.tsv`, `runs.tsv`, `issues_at_published_calls.tsv`, `published_calls_at_issues.tsv`, `kmer_databases.tsv` | what was reported without a prediction; what was counted | `tables.py` |
| `plant.sh`, `planted_*.tsv` | score 3 | `scripts/plant_errors.py`, `scripts/eval_planted.py`; copied by `tables.py` |
| `alt_settings.py`, `alt_settings_hg002.tsv`, `alt_settings_na19240.tsv` | the two candidate values | `alt_settings.py` |
| `hapmers_against_hprc.py`, `verkko_hapmers_*.tsv` | the added run | `hapmers_against_hprc.py`; copied by `tables.py` |
| `after.py`, `after_*.tsv` | descriptions added after the scored run was seen | `after.py` |
| `summary_*.json`, `run_*.json` | each run's summary; its inputs and every parameter it ran with | the pipeline; copied by `tables.py`, and `after.py` names the inputs under `DATA` |
