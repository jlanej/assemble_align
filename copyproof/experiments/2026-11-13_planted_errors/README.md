# Rung 2: errors planted in a truth genome

**Question.** Edit T2T-HG002 Q100 v1.1 so that it carries assembly errors of known class, size and
place, and run copyproof with HG002's real short reads (which come from the unedited genome). Is
each edit found, as the right class, and does anything else change (GUIDE §8 rung 2; decisions D1,
D4, D6)?

**Inputs.** `hg002v1.1.fasta.gz`; the chain between its haplotypes; its segmental-duplication,
centromere, HSat and known-issue annotation; HG002 NovaSeq PCR-free reads (36x). Parents were not
used, so the hapmer channel (S3) is silent here.

**Run.**

    scripts/plant_errors.py --fasta hg002v1.1.fasta.gz --chain hg002v1.1_to_other_haplotype.chain.gz \
        --sd SDs.bed --avoid centromeres.bed hsats.bed issues.v4.bed --assessable cannot_assess.bed --out plant1
    scripts/run_sample.py --out run35 --asm mat=plant1.fa --group-map plant1.groups.tsv --db sample=HG002.35x.cpkdb --reads R1 R2
    scripts/eval_planted.py --expected plant1.expected.bed --edits plant1.edits.tsv --run run35 --baseline TRUTH_RUN --out run35.eval

111 of 116 requested edits found a place (seed 20261113; 4 per kind, size and stratum; at least
3 Mb apart; on autosomes; outside satellites, known issues and the copy-number mask): 39 tandem
duplications, 8 extra sequences, 36 deletions, 12 haplotype swaps, 12 inversions and 4 inversions
rejoined through 100 N. `edits.tsv` lists them.

## Result

Found = an issue of the expected class covers at least half of the expected stretch (for point
signals: lies within 200 bp).

| Planted error | Expected issue | 2 kb | 5 kb | 10 kb | 50 kb | 200 kb |
|---|---|---|---|---|---|---|
| Tandem duplication, unique sequence | `EXPANSION` on the two copies | 2/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| ... its homologous stretch | `EXPANSION`, role `partner` | 0/4 | 1/4 | 4/4 | 4/4 | 4/4 |
| Tandem duplication, in a segmental duplication | `EXPANSION` | 1/4 | 4/4 | 4/4 | 4/4 | 3/3 |
| ... its homologous stretch | `EXPANSION`, role `partner` | 0/4 | 2/4 | 4/4 | 4/4 | 3/3 |
| Extra sequence identical to a stretch of one haplotype | `EXPANSION` on the extra sequence and on its twin | | | 7/8 | 8/8 | |
| ... the homologous stretch | `EXPANSION`, role `partner` | | | 4/4 | 4/4 | |
| Deletion from one haplotype, unique sequence | `COLLAPSE` on the homologue | 0/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| | `MISSING_HAPLOTYPE` on the homologue | 3/4 | 3/4 | 4/4 | 4/4 | 4/4 |
| Deletion, in a segmental duplication | `COLLAPSE` | 0/4 | 2/4 | 4/4 | 3/4 | |
| | `MISSING_HAPLOTYPE` | 2/4 | 2/4 | 3/4 | 3/4 | |
| Stretch replaced by the other haplotype's | `MISSING_HAPLOTYPE` | | | 3/4 | 4/4 | 4/4 |
| Inversion (two false adjacencies each) | `MISSING_JUNCTION` | | | 0/8 | 0/8 | 0/8 |
| Inversion rejoined through 100 N | `MISSING_JUNCTION` | | | | 0/8 | |

Of the found signals, 119 of 147 are high-confidence (`COLLAPSE` 22 of 25, `EXPANSION` 36 of 49,
partners 29 of 34, `MISSING_HAPLOTYPE` 32 of 39). The table is from the final classifier (the one
that came out of rung 4, see `../2026-11-13_hg002_assemblies/`).

**Nothing else changes.** The planted genome carries 656 issues. 237 belong to a planted edit:
at the edit, at a place its event is also anchored (`twin_of`), or at another copy of its sequence
(`other_copies`: the paralogues of an edit planted in a duplication). 398 of the other 419 are
the issues of the unedited genome (rung 3), 8 lie at the stretch homologous to an edit, and 13
are new: 7 `MISSING_SEQUENCE` and 2 `MISSING_ALLELE` calls, and 4 one- or two-window copy-number
calls of low confidence (`new_elsewhere.tsv`).

## What was learned

1. **A false adjacency is invisible to missing k-mers when the other haplotype is right.** At an
   inversion breakpoint the reads carry the true adjacency, and the release still holds it: on the
   homologue. No read k-mer is missing. What does show is the release's own false adjacency: its
   k-mers have no read. All 24 inversion breakpoints and 34 of 36 deletion joins are runs of
   unsupported assembly k-mers (`unsupported_clusters.tsv`), 23 and 29 of them with normally
   covered flanks. But rung 3 found 1,123 such runs on the unedited genome, nearly all in sequence
   that long reads confirm. So an unsupported run is reported only with a second witness: a
   missing allele from the reads (`BASE_ERROR`), or long reads that also break there
   (`FALSE_JOIN`, needs `--long-read-bam`; there is no long-read alignment to this edited genome,
   so the class is silent here). A join through N has no k-mers at all and is the business of the
   join audit (S4).
2. **An identical extra copy cannot be settled from its own k-mers.** Where the release holds a
   sequence three times and the reads show two copies, "one copy too many" and "two thirds of the
   depth" predict the same counts unless k-mers of another multiplicity lie in the stretch. A
   10-kb copy has about 80 tracked heterozygous k-mers, too few. First result: 1 of 8 found.
   Two things fixed it, both now part of the method:
   - *The other copies decide.* The k-mers of the segment are located in the release; the
     k-mers unique to those other copies lie between the shared ones and must be depleted too if
     sequencing depth is low there. Read at full depth, they refute the depth explanation.
   - *Long and deep ambiguous segments are rare in a correct genome.* On the unedited genome no
     ambiguous segment in 5.3 Gb is at least 8 kb long with a depth change of 1.5-fold or more
     (one, under the classifier of that moment); dips to half the depth that last one or two
     windows number 105. Such segments are now reported, at low confidence.
3. **An event must not outlive its evidence.** A hidden state for "one copy fewer, among k-mers
   the release holds three times or more" says nothing about a window without such k-mers; its
   likelihood there equals the concordant one, and the model stayed in it for free. Segments ran
   on past the duplicated stretch (to the end of the sequence in a simulated genome) and were then
   classified as irregular. Segments are now built from windows that test their state. This
   raised partner detection at 10 kb from 2/4 to 4/4 and cut irregular segments on the unedited
   genome from 958 to 848 (401 after the classifier revision of rung 4).
4. **Size limits.** Copy-number changes below 5 kb are mostly missed at 36x with 1-kb windows
   (the segment must gather 20 nats over a prior of 10^-6 per window); a missing haplotype is
   still found at 2 kb from its alleles (5 of 8).
5. **Boundaries matter at 10 kb.** The one extra-sequence entry that is missed (the copy it
   repeats, chr10 paternal 77.57 Mb) is a segment that runs 1.5 kb past the copy: its flank's
   two-copy k-mers mix with the copy's and the whole-copy pattern no longer fits within 4%.

## Files

- `sensitivity.tsv`: found and planted per kind, signal, stratum and size.
- `edits.tsv`: every planted edit in original and edited coordinates.
- `new_elsewhere.tsv`: the 13 issues that are neither at an edit, nor at its homologue, nor in the
  unedited genome.
