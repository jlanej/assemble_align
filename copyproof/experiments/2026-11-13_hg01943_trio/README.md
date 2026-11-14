# Rung 1: HG01943, the case that started this

**Question.** The HPRC release-2 assembly of HG01943 carries the failure this tool was designed
for (report F4): the child has a de novo 1.09-Mb tandem duplication on the maternal haplotype; the
maternal contig holds it once, and the second copy was scaffolded into the paternal chromosome
through four 100-N joins, with two false duplications at the joins. Run on the published
assembly with public reads, does copyproof flag the five defects without flooding the locus
(GUIDE §8 rung 1; regression test R2 of `../../../litmus/cases_HG01943_HG01103/`), and what does it
report on the rest of that genome?

**Inputs.** `HG01943_pat_hprc_r2_v1.0.1.fa.gz` and `HG01943_mat_hprc_r2_v1.0.1.fa.gz` (5.98 Gb);
the NYGC 30x CRAMs of the child (HG01943), mother (HG01942) and father (HG01941); HPRC's own
HiFi-to-assembly alignment, read remotely by region.

**Run.** `scripts/run_hg01943.sh DATA` (`copyproof run` with `--parent` twice and
`--long-read-bam`), then the R2 scorer on `issues.bed`.

## Result at the locus

| Target of the regression test | copyproof | Found |
|---|---|---|
| T1: the four joins that no read spans | `UNSUPPORTED_JOIN` at all four (0 linking reads, 4 to 11 reads on each flank), each high-confidence because it lies next to another issue | 4 of 4 |
| T2: the piece on the wrong haplotype | `MISPLACED_PIECE` paternal 47,396,531-47,765,590: 281 maternal hapmer sites against 0 paternal; the only misplaced piece in the genome | 1 of 1 |
| T3: the collapsed copy | `COLLAPSE` maternal 47.368-47.455 Mb (`locus`) and 47.684-47.782 Mb (`partner`), paternal 48.103-48.201 Mb (`locus`) | 2 of 2 |
| T4: the junction of the duplication | `MISSING_SEQUENCE` at maternal 46,997,978, one base from the junction end; 4 reads, k-mers in neither parent | 1 of 2 ends |
| T5: the two false duplications at the joins | `EXPANSION` paternal 47.376-47.397 and 47.745-47.786 Mb (the first pair), 48.202-48.241 and 48.524-48.541 Mb (the second) | 2 of 2 |

5 of 5 targets, 19 calls in the locus, 2 off-target (one `BASE_ERROR`, one `MISSING_SEQUENCE`):
**pass** by the rule fixed in the regression test (5 of 5 and at most 2 off-target). For
comparison, the published QC of the same assembly scores 0 to 2 of 5 (HMM-Flagger HiFi and ONT,
NucFlag).

**This is not an independent test for two of the rules.** The first run scored 3 of 5: the
junction and the second false duplication were seen but not reported. Two rules were added after
looking at them, both stated generally and both documented in `docs/method.md`:

1. *Missing sequence that is in neither parent.* The junction reads were placed within 1 bp, but
   as a weak one-sided cluster of 4 reads: the far side of the junction is a high-copy repeat,
   where reads are deliberately not recruited. Weak one-sided clusters are not reported (nine in
   ten are artefacts). But recurrent artefacts recur in the parents' reads; clusters whose k-mers
   are solid in the child and in neither parent are few (45 events in this genome), and the
   junction is the one with the highest k-mer count. These are now reported at low confidence.
2. *An ambiguous copy-number segment at an unlinked join.* The second false duplication is "one
   copy fewer among four", which a depth factor of 0.81 explains as well. It lies within 5 kb of
   a join that no long read links, and that is reported now (2 segments in this genome, both
   here).

## The rest of the genome

2,867 issue rows (`issues_by_class.tsv`). Depth per copy 14.2; 88.9% of bases assessable.

| Class | High | Low | Note |
|---|---|---|---|
| `MISPLACED_PIECE` | 1 | 0 | the known one |
| `PHASE_SWITCH` | 19 | 137 | 12.7 Mb switched in all (0.2% of the genome) |
| `UNSUPPORTED_JOIN` | 11 | 23 | 34 of the 42 joins have no linking read; 11 lie next to another issue |
| `COLLAPSE` | 14 (+8 partner, +14 unresolved) | 4 | |
| `EXPANSION` | 17 (+32, +37) | 47 | |
| `UNSUPPORTED` | 1 | 2 | |
| `FALSE_JOIN` | 5 | 3 | no read of either kind holds the k-mers, and no long read runs through |
| `BASE_ERROR` | 425 | 284 | no read of either kind holds the k-mers; long reads run through |
| `MISSING_HAPLOTYPE` | 25 (+19 twin) | 21 (+21) | |
| `MISSING_ALLELE` | 558 (+527 twin) | 198 (+189) | |
| `MISSING_JUNCTION` | 1 | 0 | |
| `MISSING_SEQUENCE` | 65 (+72 twin) | 43 (+44) | 87 rows have k-mers in neither parent |

1. **Seven more joins look like the known four.** Besides the four at the locus, seven unlinked
   joins have a copy too many right next to them (`joins.tsv`, `next_to=EXPANSION`): on
   JBHIIK010000006.1 at 22.40 Mb, JBHIIK010000076.1 at 22.89 Mb, CM089367.1 at 49.41 Mb (also a
   phase switch), JBHIIL010000001.1 at 0.83 Mb, JBHIIL010000004.1 at 6.17 Mb and JBHIIL010000072.1
   at 5.78 and 90.31 Mb. These are candidates for the same kind of error and were not examined
   further.
2. **An unlinked join alone says little.** 34 of 42 scaffold joins have no long read across them,
   which is what a gap is. Reporting each as an error (as a clip-based tool does at all 42) would
   bury the four that matter.
3. **Long reads sort the unsupported runs here too.** Of 2,184 runs of assembly k-mers that no
   short read holds, HiFi reads carry the k-mers of 1,072 and not of 659 (452 unclear). The 659
   are sequence in the release that neither read type supports. They make most of the 709
   `BASE_ERROR` rows (the others are places where short reads carry another allele) and the 8
   `FALSE_JOIN`. Of 51 copy-number losses on a haplotype's own k-mers, long reads carry 26
   (dropped or demoted) and do not carry 3.
4. **Cost.** Counting three 30x CRAMs through samtools took 29 to 35 minutes each on a busy
   machine; the scan 15 minutes on a throttled disk (93 s otherwise); recruiting reads from the
   CRAM 13 minutes; 2,235 remote long-read fetches 20 minutes. A re-run on cached stages takes 2
   minutes.

## Limits seen here

- The duplication's junction is detected only weakly, because one side of it is a high-copy
  repeat. A junction between two stretches of unique sequence is placed on both sides (rung 0).
- Three of the four classes at this locus need more than the child's short reads: the piece needs
  the parents, the joins need long reads, the junction (as it turned out) needs the parents.
  With short reads of the child alone the locus would show the collapses and the first false
  duplication.

## Files

- `score_R2.json`: the scorer's output. `issues_at_locus.tsv`: the 17 issue rows at the locus.
- `issues_by_class.tsv`, `joins.tsv` (the 34 unlinked joins with what lies next to them),
  `summary.json`.
- `card_misplaced_piece.png`, `card_collapse.png`, `card_expansion_at_join.png`: evidence cards.
