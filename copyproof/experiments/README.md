# Experiments: an index, and where each number comes from

One folder per experiment. Each holds a `README.md` (question, inputs, commands, result, what was
learned) and the small tables behind its numbers. Large files stay outside git; `DATA` in the
commands is the data root (`scripts/manifests/` lists what to download into it).

| Folder | Rung of the guide's ladder | Question |
|---|---|---|
| [`2026-11-13_engine_validation`](2026-11-13_engine_validation/README.md) | | Are the k-mer counts right, and what do they cost? |
| [`2026-11-13_hg002_q100_false_flags`](2026-11-13_hg002_q100_false_flags/README.md) | 3 | What is flagged on a truth genome, with short reads and with long reads? |
| [`2026-11-13_hg002_trio`](2026-11-13_hg002_trio/README.md) | 2 and 3, with parents | What does the hapmer channel flag on a truth genome, and does it find planted swaps? |
| [`2026-11-13_planted_errors`](2026-11-13_planted_errors/README.md) | 2 | Are errors planted in the truth genome found, and does anything else change? |
| [`2026-11-13_hg002_assemblies`](2026-11-13_hg002_assemblies/README.md) | 4 | On real assemblies of the same person, are the issues real and how many are found? |
| [`2026-11-13_hg01943_trio`](2026-11-13_hg01943_trio/README.md) | 1 | Is the case that motivated the tool flagged, and what else is in that genome? |
| [`2026-11-13_threshold_sensitivity`](2026-11-13_threshold_sensitivity/README.md) | | Do the results move with the thresholds? What did the classifier's revision change? |

Rung 0 (a synthetic trio) is the test suite, `tests/`. Rungs 5 and 6 (other samples; thresholds
frozen beforehand) were not run.

## Reading order for a reviewer

1. [`../README.md`](../README.md), "What it finds": the claims.
2. This page: where each claim's number comes from.
3. The experiment records, in the order rung 3, trio, rung 2, rung 4, rung 1: that is the order
   in which the method was shaped. The sensitivity record last.
4. [`../docs/method.md`](../docs/method.md) for the method as built and
   [`../docs/thresholds.md`](../docs/thresholds.md) for every numeric choice and its origin.
5. [`../NOTEBOOK.md`](../NOTEBOOK.md) for the order of events. Its entries were written the same
   evening, not as each experiment ran, and three of them carry counts from before later changes;
   each says so at its head. Where the notebook and a record differ, the record is regenerated
   by a script and is the reference.

The code was written and changed in one day and arrives as one large commit, so the earlier
versions of the method are not in the history. Two things stand in for that: the classifier as
it was before its revision is kept in `scripts/dev/classifier_before_revision.py`, and the
sensitivity record runs both on the same segments; the hapmer audit can be run with and without
its two filters (`2026-11-13_hg002_trio/tables.py`).

## Claims and their sources

Every number in `../README.md`, with the table it is read from and what writes that table.

| Claim | Table | Written by |
|---|---|---|
| High-confidence copy-number issues confirmed by the truth genome; assessable discordances found, per assembly and summed | `2026-11-13_hg002_assemblies/scoreboard.tsv`; the sums as quoted are written out in `scoreboard.md` | `scripts/eval_truth.py` per assembly, then `scripts/scoreboard.py` |
| Found and confirmed by size (surplus copies of 5 kb and more, collapses, under 5 kb) | `2026-11-13_hg002_assemblies/scoreboard.by_size.tsv` | the same |
| The unconfirmed high-confidence issues, 1 to 10 per assembly, and where they lie | `2026-11-13_hg002_assemblies/unconfirmed_high.tsv` | `2026-11-13_hg002_assemblies/unconfirmed.py` |
| Strong missing-allele calls true 96% of the time, one-sided missing sequence 46% | `2026-11-13_hg002_assemblies/hprc_year2.missing_clusters.tsv` | `2026-11-13_hg002_assemblies/details.py` (the run needs `WITH_READS=1`) |
| Of the collapses missed in one assembly, how many show no departure in the reads | `2026-11-13_hg002_assemblies/hprc_year2.collapses.tsv` | `details.py` |
| Whole-copy calls under 5 windows are mostly not confirmed | `2026-11-13_hg002_assemblies/hprc_year2.whole_copy_by_length.tsv` | `details.py` |
| Planted errors found, by kind and size | `2026-11-13_planted_errors/sensitivity.tsv` | `scripts/eval_planted.py` |
| Planted false adjacencies show as unread assembly k-mers | `2026-11-13_planted_errors/breakpoints.tsv` | `scripts/eval_planted.py` |
| Issues on the planted genome that are neither at an edit nor on the unedited genome | `2026-11-13_planted_errors/new_elsewhere.tsv` (the totals are printed) | `scripts/eval_planted.py` |
| Issues on the truth genome by class and confidence (about 25 high-confidence per haplotype) | `2026-11-13_hg002_q100_false_flags/issues_by_class.tsv` | `2026-11-13_hg002_q100_false_flags/analyse.py` |
| 0.11% of 20,237 duplication boundaries flagged | `2026-11-13_hg002_q100_false_flags/litmus.txt` | `scripts/sd_boundary_litmus.py` (the denominator was fixed in NOTEBOOK entry g before the run) |
| No ambiguous segment on the truth genome is long and deep; 105 are short dips | `2026-11-13_hg002_q100_false_flags/ambiguous_by_length.tsv` | `analyse.py` |
| Long reads hold the k-mers of 957 of 1,123 unsupported runs; 115 of the other 119 are known issues, 27 changed in v1.2 | `2026-11-13_hg002_q100_false_flags/unsupported_runs_hifi.tsv` | a run with `--long-read-bam`, then `annotate_runs.py` |
| What the reads hold instead at those 119 (106 of 110 a different length) and at the four on no list | `not_carried_reads_allele.tsv`, `candidates_reads_allele.tsv` in the same folder | `reads_allele.py` |
| Long reads hold the unread k-mers at 10 of 13 short copy-number losses | `2026-11-13_hg002_q100_false_flags/hifi_at_short_losses.tsv` | `hifi_at_short_losses.py` |
| Hapmer channel on the truth genome: 0 of 18 pieces, 143 then 24 switch segments | `2026-11-13_hg002_trio/truth_filters.tsv`, `truth_switches.tsv` | `2026-11-13_hg002_trio/tables.py` |
| Planted swaps found by the hapmer channel (7 of 12) | `2026-11-13_hg002_trio/planted_swaps.tsv` | `tables.py`; the 11 of 12 by the read channel are in the planted `sensitivity.tsv` |
| HG01943: 5 of 5 targets, 2 off-target | `2026-11-13_hg01943_trio/score_R2.json` | `../../litmus/cases_HG01943_HG01103/score_flagged_cases.py` on the run's `issues.bed` |
| HG01943: one misplaced piece; 34 unlinked joins, 11 beside another issue; none of the 7 outside the locus repeats sequence from another join | `2026-11-13_hg01943_trio/issues_by_class.tsv`, `joins.tsv`, `joins_evidence.tsv` | `2026-11-13_hg01943_trio/tables.py` |
| HG01943: missing sequence in neither parent is 46 events, the known junction first | `2026-11-13_hg01943_trio/missing_in_neither_parent.tsv` | `tables.py` |
| The engine's counts against jellyfish | `2026-11-13_engine_validation/compare_jellyfish.tsv` | `compare_jellyfish.py` |
| The classifier before and after its revision; one threshold moved at a time | `2026-11-13_threshold_sensitivity/sensitivity.tsv` | `scripts/dev/cache_segments.py`, then `scripts/dev/sensitivity.py` |
| Costs (minutes, GB) | stated in `../README.md` and in the engine record | timings printed by the engine and in each run's `summary.json`; not collected into a table |

## What cannot be regenerated from this repository

- The intermediate states of the day. The notebook's entry j quotes results of pipeline runs made
  between changes (for example "collapses found 8 of 64, then 17"); those runs were overwritten.
  The sensitivity record gives the before and after of the classifier on cached segments, which
  agrees with them (8 and 17 of 64) but leaves out two pipeline steps.
- The first exploratory runs at half depth, which the notebook's entries d to f describe.
- Anything at all without the downloads: about 220 GB of reads, assemblies and alignments
  (`scripts/manifests/`). The k-mer databases (32 GB) and the run directories (about 315 GB as
  they stand, nearly all of it positional tracks of 9 to 18 GB per run) are rebuilt from those.
