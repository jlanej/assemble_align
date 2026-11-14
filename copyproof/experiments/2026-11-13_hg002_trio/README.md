# The hapmer channel on a truth genome: HG002 with its parents

**Question.** The parental-hapmer channel (S3) had been tested on a simulated trio and on one real
case (HG01943). What does it flag on a genome whose phasing is right, and does it find haplotype
swaps of known size (GUIDE §8 rungs 2 and 3; decision D7)? And what do the parents' reads say
about the read sequence that the truth genome lacks?

**Inputs.** T2T-HG002 Q100 v1.1 and the same genome with 111 planted errors (12 of them stretches
of 10, 50 or 200 kb replaced by the other haplotype's sequence); NovaSeq PCR-free reads of HG002,
of the mother HG004 and of the father HG003 (about 36x each).

**Run.** `scripts/run_hg002_trio.sh DATA` (counts the parents' k-mers once, 5 minutes each; then
both genomes with three read sets), then `tables.py` here, which repeats the audit with and
without its two filters on the same windows.

## Result

**On the truth genome every phase switch is a false flag.** All 46 chromosomes are assigned to
the right parent, with 98.9% to 99.97% of their hapmer sites agreeing, and none of the 18
gap-delimited pieces is called misplaced. Switch segments:

| Audit | Switch segments | ... with 10 sites or more | Mb |
|---|---|---|---|
| as first written | 143 | 28 | 10.9 |
| only where both parents are readable | 32 | 8 | 5.6 |
| only switches with sites at a tenth of the usual density | 125 | 25 | 1.4 |
| both (as shipped) | 24 | 8 | 0.4 |

**Planted swaps found as a phase switch** (at least half of the stretch covered): 1 of 4 at 10 kb,
2 of 4 at 50 kb, 4 of 4 at 200 kb, the same with and without the filters. The missing-haplotype
call of the read channel, which needs no parents, finds 11 of the same 12.

## What was learned

1. **Most false switches were places a parent's reads do not cover.** A hapmer says "this
   parent's reads hold the k-mer and the other's do not". Where the other parent's reads are
   thin, the second half means nothing. The largest false switches lay at the immunoglobulin
   loci (chr14 107.7 to 108.3 Mb, chr2 90.2 Mb, chr22 25.0 to 25.4 Mb, all on the maternal
   haplotype, up to 151 "paternal" sites against 7): the mother's sample is a B-cell line, which
   has deleted part of the haplotype she passed on, and the father happens to carry the same
   alleles. Others were ordinary dropouts (the GC-rich stretch at chr2 1.53 Mb, acrocentric
   short arms). The audit now sets aside windows where a parent's depth on the two-copy k-mers
   is outside 0.7 to 1.5 of its usual depth: 1.7% of windows, and 111 of the 143 segments.
2. **A few stray sites over a megabase are not a switch.** Seven segments of 0.7 to 1.4 Mb held
   5 to 14 sites each: stretches where the parents share haplotypes, so that almost no site says
   which parent the sequence came from, and the half percent of sites that are wrong decide. A
   switch must now carry the other parent's sites at a tenth or more of the sequence's own
   density. Any value from 0.1 to 0.2 gives the same 24 segments; at 0.2 one planted 200-kb
   swap in a stretch poor in sites (15 sites, ratio 0.20) is lost.
3. **What is left.** 24 segments in 0.4 Mb; 8 have 10 sites or more, and 4 of those are still
   the immunoglobulin heavy-chain locus, where the mother's depth stays above the bound
   (`truth_switches.tsv`). So a phase switch at an immunoglobulin locus of a cell-line trio says
   something about the parents' cells, not the assembly.
4. **Short switches are the read channel's business.** A 10-kb swap holds about four hapmer
   sites, at the edge of what the channel accepts. The missing-haplotype call sees the same swap
   as a run of alleles the reads carry and the release lacks, without parents. S3 earns its
   place with pieces of hundreds of kb (the HG01943 piece: 281 sites against 0).
5. **A parent's reads do not vouch for missing sequence.** A rule written before this run made a
   strong one-sided missing-sequence call high-confidence when a parent's reads held its k-mers.
   On the truth genome a parent holds the k-mers of 208 of 211 strong clusters: recurrent
   artefacts recur in every library. The rule is gone; such calls are low-confidence. The
   opposite use of the parents holds up: missing k-mers that are solid in the child and in
   neither parent are 9 placements (5 events) on the truth genome, against 46 events in HG01943
   where the known de novo junction is first among them.
6. **None of this tests the truth genome's phasing.** Its haplotypes were assigned with parental
   data (they carry the parents' names), so 46 chromosomes on the right parent is agreement
   between two readings of the same kind of evidence. What this run measures is the channel's
   false flags. The audit is an independent test of phasing only on an assembly phased without
   the parents; on a trio-phased one it finds errors made after phasing, such as the HG01943
   piece.

## Decisions

- **D7:** S3 needs both parents, and needs them readable; it is specific for misplaced pieces (0
  of 18 on the truth genome, 1 of 1 in HG01943) and weak for switches under 50 kb.
- **D6:** the two filters and the removal of the inheritance rule. As with the other channels,
  the thresholds were set on the data they are reported on.

## Files

- `tables.py`: repeats the audit under the four settings on both genomes and prints the tables.
- `truth_filters.tsv`, `planted_swaps.tsv`: the two tables above.
- `truth_switches.tsv`: the 24 segments left on the truth genome.
