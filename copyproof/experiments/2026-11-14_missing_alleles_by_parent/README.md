# Alleles in the reads and not in the release, by the parents' reads

**Question.** The held-out test found that half of NA19240's strong missing-allele calls have
their k-mers in neither parent's reads, against 7% on HG002, and did not look further
(`../2026-11-14_na19240_heldout/README.md`, "Missing alleles"). Is that NA19240 alone? Did the
parents' reads fail to cover sequence that she inherited? And what do runs that already exist say
about what these calls are?

This record tabulates finished runs. Nothing was run again, no rule was made and no code under
`src/` changed. It was made after the held-out test and from its runs, so nothing here is held out.

**Inputs.** Six runs, each with the child's short reads and both parents' (`runs.tsv`):

| Run | Genome audited | Reads |
|---|---|---|
| `HG002_hprc_year2` | HG002, the assembly `hprc_year2` | NovaSeq: HG002, HG004, HG003 |
| `HG002_truth_genome` | HG002, T2T Q100 v1.1 | the same |
| `HG01943_release` | HG01943, HPRC release 2 | NYGC: HG01943, HG01942, HG01941 |
| `NA19240_release` | NA19240, HPRC release 2 | NYGC: NA19240, NA19238, NA19239 |
| `HG002_truth_genome_planted` | the truth genome with 111 planted errors | as the truth genome |
| `NA19240_release_planted` | NA19240's release with 115 planted errors | as the release |

The depth of one copy is 13.4 to 14.2 and the channel's threshold is 7 reads in every run. The
runs of NA19240 and of `hprc_year2` are those of the held-out test (commit `1465496`). The other
three were made on the evening of 13 November, before that commit, and were not made again.
`src/copyproof/missing.py`, which writes the cluster table read here, has a single commit, older
than all six runs. The step from clusters to issues changed once that evening, for one-sided
sequence only, and both older issue lists already follow the later rule (no one-sided call is of
high confidence in `issues_by_parent.tsv`).

**Run.**

    .venv/bin/python experiments/2026-11-14_missing_alleles_by_parent/tables.py DATA experiments/2026-11-14_missing_alleles_by_parent

It reads the run directories and the side files of the read databases, writes nothing under
`DATA`, and takes 4 seconds.

**Terms.** A *placement* is one row of a run's `missing_clusters.tsv`: reads that agree on one
place of one sequence (the held-out record says "call"). An allele whose flanks exist on both
haplotypes is placed on both, so one *event* is usually two placements (`twin_group`). *Strong*
is what makes a missing allele a high-confidence issue in `calls.py`: 6 reads or more and 5
tracked k-mers or more. *In a parent's reads*: half or more of the placement's tracked k-mers are
in the mother's read database, or half or more in the father's. A read database keeps a k-mer
from its second sighting, so "in" means seen at least twice. The *read count* of a placement is
the median count of its k-mers in the child's reads.

## Result

**Every run has such placements, and the two HPRC releases have hundreds**
(`clusters_by_parent.tsv`). Strong missing-allele placements:

| Run | All | In a parent's reads | In neither | Events: in a parent, in neither | Median read count: in a parent, in neither |
|---|---|---|---|---|---|
| HG002, `hprc_year2` | 996 | 925 | 71 (7%) | 469, 38 | 13, 8.5 |
| HG002, truth genome | 74 | 40 | 34 (46%) | 26, 18 | 7, 8 |
| HG01943, release | 1,701 | 1,156 | 545 (32%) | 731, 276 | 14.5, 8 |
| NA19240, release | 1,773 | 851 | 922 (52%) | 467, 460 | 14, 8 |

HG01943, the sample of rung 1, shows what NA19240 shows. Among the weak placements, 17, 8, 197
and 298 are in neither parent's reads.

**The split is all or nothing.** Of the strong placements in neither parent's reads, 71 of 71,
34 of 34, 545 of 545 and 911 of 922 have not one k-mer in either parent's database. Of those in a
parent's reads, 921 of 925, 34 of 40, 1,149 of 1,156 and 818 of 851 have every k-mer in one. A
parental share between 0.25 and 0.75, where the cut at one half could matter, occurs in 0, 0, 4
and 17 strong placements (`checks.tsv`).

**The parents' reads did not miss them.** Two measurements say so.

- *The reads* (`read_sets.tsv`). All nine read databases keep a k-mer from its second sighting.
  The parents of the two NYGC trios have the same amount of reads as HG002's parents or more:
  141.0 and 108.6 Gb for NA19240's, 116.3 and 119.2 Gb for HG01943's, 108.8 and 108.5 Gb for
  HG002's. The peak of each parent's count histogram lies at 26 to 34, so one copy gives 13 to
  17 reads. At that depth a k-mer present once in a parent's genome is seen fewer than twice with
  a chance of less than 4 in 100,000 (Poisson; the table gives it for each database).
- *A control* (`planted_control.tsv`). Errors planted into a genome take alleles out of it: among
  the 111 and 115 edits are 36 and 39 deletions from one haplotype and 12 stretches replaced by
  the other haplotype's sequence. The planted runs hold 2,019 and 2,658 strong placements that the
  unedited runs do not (1,422 and 1,846 events). All 4,677 are in a parent's reads. None is in
  neither, which bounds that share at 0.064% (one-sided, 95%). Every strong placement of the
  unedited runs is found again in the planted ones, so the 34 and the 922 are the same before and
  after the edits.

An allele that the child inherited and the release lacks is found in a parent's reads. The
placements in neither parent's reads are not inherited sequence that went unread.

**Their read counts are below one copy** (`alleles_by_read_count.tsv`,
`read_counts_against_planted.tsv`). The placements made by planted errors have a median read count
of 14 and a lower quartile of 11. Below that quartile lie 798 of NA19240's 922 placements in
neither parent's reads (87%), 54 of 71 on `hprc_year2` (76%) and 24 of 34 on the truth genome
(71%). HG01943 has no planted run; against its own placements in a parent's reads (lower quartile
11.5) the number is 509 of 545 (93%). A read count of 13 or more is found in 11 of 71, 4 of 34,
10 of 545 and 24 of 922. The counts begin at the threshold of 7, and 20, 11, 155 and 245
placements lie under 8: the threshold cuts into this kind of placement, and what is counted here
is its upper part.

**They look like single changed bases, read from both strands and scattered**
(`clusters_by_parent.tsv`).

- The cluster table's `size` is 0 at 71 of 71, 34 of 34, 543 of 545 and 919 of 922: between the
  flanks the reads' length differs from the release's by 30 bases or fewer. Corrected on 14
  November 2026: this bullet first said "the same length as the release". `size` is 0 whenever one
  chain of anchors spans the missing k-mers, and a chain takes in a length difference of up to 30
  bases (`diag_tol` in `src/copyproof/missing.py`). The record that followed derives the change
  from the reads: one changed base at 70 of 71, 32 of 34, 506 of 545 and 826 of 922 placements,
  and another length at 0, 2, 18 and 39
  (`../2026-11-14_missing_alleles_long_reads/changed_bases.tsv`).
- The median is 7 tracked k-mers in every run. One changed base makes 31 k-mers that the release
  lacks, and a quarter of all k-mers are tracked.
- The rarer strand carries a median of 38% to 41% of a placement's reads. For placements in a
  parent's reads it is 33% to 40%.
- Within 10 kb of another placement of the kind lie 0 of 71, 0 of 34, 16 of 545 and 20 of 922. Of
  the placements in a parent's reads, 290 of 925, 10 of 40, 648 of 1,156 and 526 of 851 have such
  a neighbour. No megabase holds more than 5 of the first kind; one holds 258 of the second.

**Most are real sequence of the person's cells** (`alleles_held.tsv`,
`alleles_assemblies_agree.tsv`).

- HG002, `hprc_year2`: the truth genome holds the k-mers of 51 of the 71 (72%), against 906 of
  925 for placements in a parent's reads (98%). Each of four other assemblies of HG002 holds 50 to
  58 of the 71.
- NA19240: one or both HGSVC3 assemblies hold 621 of the 922 (67%), against 774 of 851 (91%). The
  two assemblies were built from the same reads. They disagree on 259 of the 922 (28%: one holds
  the allele and the other does not) and on 61 of the 851 (7%).
- HG01943: no other assembly of this person is at hand.
- The 34 placements on the truth genome (18 events, `truth_genome_alleles.tsv`) are read sequence
  that the truth genome itself lacks.

**Where a truth genome exists, the read count separates the placements about as well as the
parents do, and the two overlap** (`alleles_held.tsv`). On `hprc_year2`, placements whose k-mers
the truth genome holds:

| Read count | In a parent's reads | In neither | Together |
|---|---|---|---|
| under 9 | 75 of 88 | 21 of 38 | 96 of 126 (76%) |
| 9 to 12.5 | 303 of 309 | 19 of 22 | 322 of 331 (97%) |
| 13 or more | 528 of 528 | 11 of 11 | 539 of 539 |
| all | 906 of 925 (98%) | 51 of 71 (72%) | 957 of 996 (96%) |

Of the 39 strong placements that the truth genome does not hold, 30 have a read count under 9
and 20 are in neither parent's reads. On NA19240 the two other assemblies order the groups
another way: under 9 reads they hold 357 of 489 placements in neither parent's reads (73%) and 39
of 70 in a parent's (56%); from 13 reads, 16 of 24 and 500 of 520.

**In the issue lists these placements are much of what a reader sees** (`issues_by_parent.tsv`).
Three or more alleles less than 10 kb apart are reported once, as a missing haplotype, and all
but 9 of the 321 such issues in the four runs are in a parent's reads. What is left as single
high-confidence `MISSING_ALLELE` issues leans the other way:

| Run | High-confidence `MISSING_ALLELE` issues | Of them in neither parent's reads | Issues in the whole list | Of them from this channel and in neither parent's reads |
|---|---|---|---|---|
| HG002, `hprc_year2` | 734 | 71 (10%) | 1,442 | 102 (7%) |
| HG002, truth genome | 62 | 34 (55%) | 434 | 51 (12%) |
| HG01943, release | 1,085 | 539 (50%) | 2,737 | 825 (30%) |
| NA19240, release | 1,258 | 906 (72%) | 2,438 | 1,281 (53%) |

**The k-mers show the same, less sharply** (`kmers_by_parent.tsv`,
`kmers_by_count_and_parent.tsv`). Of the read k-mers that the channel tries to place, 1.0% are
in neither parent on both HG002 genomes, 2.0% on HG01943 and 2.6% on NA19240. By their count in
the child's reads they thin out steadily, roughly halving with each read from 7 to 10 (NA19240:
12,837, 6,323, 3,736 and 2,283 tracked k-mers; `hprc_year2`: 8,528, 3,129, 1,331 and 735). There
is no peak at the depth of one copy. Under the threshold they are common: at 3 reads, between 44%
and 48% of the missing k-mers are in neither parent.

## What was learned

1. **It is not NA19240 alone.** In the two HPRC releases a third and a half of the strong
   missing-allele placements are in neither parent's reads; in HG002's assembly of the same
   pipeline, 7%. The two releases have reads of one kind (NYGC 30x) and HG002's are another set,
   and the cell lines differ as well. These runs cannot tell the two apart.
2. **It is not a gap in the parents' reads.** The parents' reads are as deep, the databases are
   built alike, and all 4,677 placements made by planted errors are in a parent's reads.
3. **The placements behave like sequence that is in part of the child's cells and in neither
   parent**: single bases, on both strands, scattered, below one copy, mostly present in other
   assemblies of the same person, and four times as often a point on which two assemblies from
   the same reads disagree. Changes that arose in a cell line would look like this. No table
   here shows that they are that. The steady fall of the k-mer counts is what such changes give
   and also what read errors that recur give.
4. **For such a placement the statement "the reads carry an allele that the release lacks" is
   mostly true, and it does not follow that the assembler erred.** If the change is in part of
   the cells, a release that holds the other allele may be right for the DNA it was built from.
   The issue list words and ranks these placements like inherited alleles that went missing.
5. **The 96% in the main README is HG002's figure** (957 of 996), from a sample where 7% of the
   strong placements are of this kind. For the kind itself the truth genome holds 51 of 71.
6. **Read count and parents carry overlapping information.** On HG002 either would have marked
   most of the 39 placements that the truth genome does not hold. On NA19240 the stand-in orders
   the groups another way. Choosing between them needs a truth genome for a sample with many
   such placements, and no sample here has both.

## What this does not show

- What the placements are. The changed bases themselves, the long reads each release was built
  from, or a second short-read library of the same person would say more. None was looked at.
- Whether a release's own DNA carries them, and so whether any of them is an error of the
  assembly.
- Anything about HG01943's 545 beyond their number and shape.
- Two independent results on HG002: both of its rows use the same reads.
- A rule. No threshold was derived. HG002, HG01943 and NA19240 have now all been looked at for
  this question, so a rule made from these tables would need another person to test it.
- That the three runs of 13 November would come out the same today.

## Elsewhere in the repository

Two sentences were made exact with this record, and nothing else was changed:

- `../../README.md` said that strong missing-allele calls are true 96% of the time. It now names
  HG002 and gives the shares above.
- `../../docs/method.md` said that missing k-mers which are solid in the child and in neither
  parent are few, with 46 events in HG01943 and 5 on the truth genome. Those are counts of
  one-sided events. Strong alleles of that kind number 276 and 18 events in the same two runs,
  and the sentence now says so.

## Files

- `tables.py`: writes every table below from the run directories.
- `runs.tsv`, `read_sets.tsv`: the six runs and the nine read databases.
- `clusters_by_parent.tsv`: placements by kind, strength and the parents' reads, with read
  counts, shape and spacing. One-sided sequence and junctions are in it as well.
- `alleles_by_read_count.tsv`, `read_counts_against_planted.tsv`: strong missing-allele placements
  by read count, beside those that planted errors made.
- `planted_control.tsv`: the control.
- `alleles_held.tsv`, `alleles_assemblies_agree.tsv`: what the truth genome and other assemblies
  hold.
- `issues_by_parent.tsv`: the issue lists.
- `kmers_by_parent.tsv`, `kmers_by_count_and_parent.tsv`: the k-mers before any read is placed.
- `truth_genome_alleles.tsv`: the 34 placements on the truth genome, with their positions.
- `checks.tsv`: checks of the definitions (the mother is the channel's first parent; no event
  lies on both sides of the parental cut, except one on NA19240; the k-mer table agrees with each
  run's own tally).
