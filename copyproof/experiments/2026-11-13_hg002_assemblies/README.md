# Rung 4: real HG002 assemblies, scored by the truth genome's k-mers

**Question.** On assemblies with real errors, are copyproof's copy-number issues real, how many of
the real discordances does it find, and what do short-read k-mers say about consensus quality
(GUIDE §8 rung 4; decisions D4, D6)?

**Inputs.** Twelve HG002 assemblies of different years, methods and data, three of them earlier
drafts of the truth genome itself (`assemblies.tsv`; all from the public Q100 benchmark folder);
HG002's NovaSeq PCR-free reads (36x), counted once; the k-mers of T2T-HG002 Q100 v1.1, counted
exactly. The read k-mer database built for rung 3 serves every assembly. Two of the assemblies
(`verkko_2023`, `hprc_curated`) carry no haplotype in their sequence names and are run as one
group. The benchmark folder lists one more name, `hifi_q28_trio_hic`; its file is the same as
that of `verkko_asm_thic`, byte for byte, and it is counted once.

**Truth without alignment.** Each assembly is scanned with the truth genome's k-mer counts as an
extra read set, so every tracked k-mer of the assembly carries the copies the truth holds
(`copyproof/truth.py`). A *discordance* is a run of windows in which at least half of the k-mers
have more (gain) or fewer (loss) copies in the truth than in the assembly. An issue is
*confirmed* when a discordance of the same direction covers at least half of it; a discordance is
*found* when issues cover at least half of it. No aligner is involved, so the labels do not
inherit an aligner's choices in repeats. The label is about copy number of k-mers: the homologous
stretch of a collapsed copy is a true discordance too.

**Run.** `scripts/run_hg002_assemblies.sh DATA assemblies.tsv` (copy number and k-mer statistics;
`WITH_READS=1` adds the missing-sequence channel), `scripts/eval_truth.py`, `scripts/scoreboard.py`.

## Result

![Scoreboard](scoreboard.png)

| Assembly | What it is | Size (Gb), sequences | Discordant with the truth (Mb) | Assessable discordances found | High-confidence issues confirmed | Low-confidence issues confirmed | QV by truth k-mers / by read k-mers |
|---|---|---|---|---|---|---|---|
| `verkko_2023` | Verkko HiFi + ONT, trio (2023); names carry no haplotype | 6.08, 3,258 | 149.6 | 3,802 of 4,937 (77%) | 3,380 of 3,390 (100%) | 429 of 551 (78%) | 44.3 / 42.03 |
| `hprc_year1` | HPRC year-1 assembly (hifiasm trio, 2021 data) | 6.02, 1,055 | 59.7 | 290 of 681 (43%) | 254 of 260 (98%) | 50 of 88 (57%) | 52.4 / 49.72 |
| `hprc_curated` | HPRC curated HG002 (accession names; one group) | 6.02, 869 | 59.5 | 273 of 667 (41%) | 233 of 242 (96%) | 60 of 179 (34%) | 52.6 / 49.82 |
| `hifiasm_2021` | hifiasm 0.12, HiFi only (2021) | 5.97, 1,030 | 53.4 | 350 of 838 (42%) | 225 of 234 (96%) | 85 of 135 (63%) | 50.8 / 48.43 |
| `verkko_asm_thic` | Verkko, trio + Hi-C | 5.97, 90 | 23.6 | 89 of 158 (56%) | 77 of 81 (95%) | 5 of 52 (10%) | 55.8 / 50.07 |
| `hifiasm_ontonly_2025` | hifiasm, ONT only (2025) | 5.98, 163 | 18.8 | 76 of 159 (48%) | 71 of 80 (89%) | 4 of 57 (7%) | 45.7 / 43.89 |
| `hprc_year2_v2_thic` | HPRC year-2 v2, trio + Hi-C | 5.98, 209 | 17.9 | 73 of 155 (47%) | 69 of 70 (99%) | 10 of 45 (22%) | 52.7 / 48.79 |
| `hprc_year2` | HPRC year-2 pipeline, hifiasm trio, polished | 5.98, 302 | 16.5 | 151 of 233 (65%) | 114 of 119 (96%) | 31 of 74 (42%) | 51.4 / 49.21 |
| `v2_trio` | Verkko v2, trio | 6.00, 678 | 15.7 | 71 of 178 (40%) | 69 of 74 (93%) | 7 of 57 (12%) | 50.6 / 49.54 |
| `hg002v0.7` | Q100 draft v0.7 | 6.00, 48 | 2.2 | 2 of 12 (17%) | 4 of 7 (57%) | 2 of 38 (5%) | 56.9 / 51.18 |
| `hg002v0.9` | Q100 draft v0.9 | 6.00, 48 | 2.1 | 0 of 5 (0%) | 3 of 6 (50%) | 0 of 32 (0%) | 58.7 / 51.78 |
| `hg002v1.0.1` | Q100 v1.0.1 | 6.00, 48 | 0.2 | 0 of 2 (0%) | 0 of 3 (0%) | 0 of 43 (0%) | 59.3 / 51.95 |

Summed over the 12 assemblies, by size of the segment (`scoreboard.by_size.tsv`):

| Direction | Size | Assessable discordances | Found | Issues | Confirmed | High-confidence issues | Confirmed |
|---|---|---|---|---|---|---|---|
| The truth holds more (collapse) | under 5 kb | 435 | 39 (9%) | 26 | 19 (73%) | 0 |  |
|  | 5 to 20 kb | 534 | 130 (24%) | 153 | 139 (91%) | 138 | 127 (92%) |
|  | 20 to 100 kb | 422 | 213 (50%) | 225 | 210 (93%) | 213 | 203 (95%) |
|  | over 100 kb | 158 | 85 (54%) | 148 | 134 (91%) | 143 | 131 (92%) |
| The truth holds fewer (surplus copy) | under 5 kb | 363 | 65 (18%) | 656 | 26 (4%) | 4 | 4 (100%) |
|  | 5 to 20 kb | 3,505 | 2,771 (79%) | 2,544 | 2,517 (99%) | 2,132 | 2,120 (99%) |
|  | 20 to 100 kb | 2,166 | 1,563 (72%) | 1,815 | 1,804 (99%) | 1,616 | 1,606 (99%) |
|  | over 100 kb | 442 | 311 (70%) | 350 | 333 (95%) | 320 | 308 (96%) |

## What was learned

1. **High-confidence copy-number issues are real, down to a floor of a few per genome.** 4,499 of
   4,566 are confirmed (99%); without the one assembly that contributes most of them, 1,119 of
   1,176 (95%). Among the nine assemblies with 20 such issues or more, 89% to 100% each. The
   ones not confirmed number 1 to 10 per assembly (67 in all, `unconfirmed_high.tsv`) whether
   the confirmed ones number 0 or 3,380, so a share says little about an assembly with few
   errors. The three drafts of the truth genome show that floor alone: 3 unconfirmed issues
   each, against 4, 3 and 0 confirmed. Two of the three are the same in every draft: a pair of
   calls 57 to 72 kb into chromosome 1 on both haplotypes, where the reads hold four more
   copies than the genome. The truth genome's own run makes the same two calls (rung 3). Of
   the 67, 43 lie where the truth shows no discordance at all and 21 where it shows one over
   20% to just under 50% of the issue.
2. **Surplus copies are found better than collapses.** For surplus copies of 5 kb and more, 70% to
   79% of the assessable discordances are found. For collapses it is 24% to 54%. A surplus copy
   of sequence the assembly holds three times lowers the depth of its k-mers by a third; a
   collapse among six copies raises it by a sixth, which is the size of the regional wander of
   depth. Most missed collapses are of that kind: in duplicated sequence where the assembly holds
   4 to 8 copies. Of the 64 assessable collapses of `hprc_year2`, 17 are found; of the other 47,
   20 show no departure in the reads at all (k-mers the assembly holds 6 times are read at 6.1
   copies where the truth holds 7: there the short reads agree with the assembly, not with the
   truth genome), 20 are seen as a change of depth by one factor, and 7 are covered in part or left ambiguous
   (`hprc_year2.collapses.tsv`, from `details.py`).
3. **Under 5 kb the channel is weak.** 9% to 18% of the discordances are found, and 4% of the
   surplus-copy issues are confirmed; nearly all of these issues are of low confidence by rule
   (fewer than 5 windows). This limit comes from 1-kb windows at 36x and a prior of 10^-6 per
   window.
4. **The classifier was revised here.** With the classifier of rungs 2 and 3, `hprc_year2` gave 8
   of 64 collapses found. The missed ones were long stretches where k-mers the assembly holds 4
   and 6 times were read at 4.9 and 7.3 copies while the 1- and 2-copy k-mers were exactly
   concordant: a whole-copy pattern, called ambiguous because each k-mer class was allowed the
   full regional wander of depth (7.7%) on its own. That wander is shared by all k-mers of a
   stretch. Giving every hypothesis one regional factor for all classes, and each class 4%
   around it, found 17 of 64 with 22 of 24 high-confidence collapses confirmed, kept the truth
   genome at 20 whole-copy segments (14 before), and lowered the planted-error results nowhere
   after a second change: classes that hold under 10% of a segment's k-mers may not veto a fit,
   because they are mixtures (a k-mer lying twice in a surplus copy loses two copies). Those
   counts are from pipeline runs between the changes and cannot be run again; the classifier of
   before is kept, and `../2026-11-13_threshold_sensitivity/` runs both on the same segments
   (collapses 8 and 17 of 64 there too, confirmed high-confidence issues 89 and 114).
5. **Short whole-copy calls are mostly not errors.** Of 18 whole-copy surplus calls on fewer than
   5 windows in `hprc_year2`, the truth confirms one (`hprc_year2.whole_copy_by_length.tsv`). They are the same family as 13 calls on
   the truth genome itself: 1 to 3 kb where one haplotype's own k-mers have no short reads. Long
   reads showed that the sequence is real in 10 of those 13 (rung 3). Such calls are now of low
   confidence, and dropped when long reads carry the k-mers.
6. **Depth anomalies are rightly left out.** The twelve runs hold 32,740 segments that fit a depth
   factor or nothing regular. The truth confirms 1,302 (4%); without `verkko_2023`, 440 of 29,305
   (2%). Ambiguous segments without support: 1,337, of which 171 confirmed (13%).
7. **One release holds 150 Mb more than a diploid genome.** `verkko_2023` (3,258 sequences) has
   149.6 Mb of sequence discordant with the truth, nearly all surplus copies; copyproof raises
   3,390 high-confidence issues there and 3,380 are confirmed. The other ordinary assemblies carry
   16 to 60 Mb of discordant sequence, the three drafts 2.2, 2.1 and 0.2 Mb.
8. **Short-read k-mers understate consensus quality, by a floor.** The QV implied by assembly
   k-mers without any read is 2 to 7 units below the QV by the truth's k-mers, and the gap grows
   with quality: 1.1 to 2.8 units in the seven assemblies below QV 52.7, 3.9 to 5.7 in the three
   up to QV 57, and 6.9 and 7.3 in the two latest drafts. By read k-mers the three drafts reach
   51.2, 51.8 and 52.0 while the truth's k-mers give them 56.9, 58.7 and 59.3: 2.4 units between
   the drafts shrink to 0.8. The reason is the share of unread k-mers that are in fact correct
   sequence. In the two weakest assemblies most unread k-mers are errors (`verkko_2023` 90%
   wrong by the truth, `hifiasm_ontonly_2025` 86%); in the seven ordinary ones 40% to 62%; in
   the drafts 24%, 13% and 9%. The reads find 32% to 92% of the k-mers that are wrong by the
   truth.
9. **Missing alleles are real; one-sided missing sequence mostly is not.** For `hprc_year2` the
   missing-sequence channel was run with the reads (`hprc_year2.missing_clusters.tsv`). A cluster
   is true if the truth genome holds most of its missing k-mers. `MISSING_ALLELE`: 957 of 996
   strong clusters are true (96%), 175 of 273 weak ones. `MISSING_SEQUENCE` (reads anchored on one side only): 53 of 115 strong ones
   (46%), 92 of 1,010 weak ones (9%). The planted genome gives the same picture (96%, 64%, 51%,
   10%). Hence weak one-sided clusters are not issues, and strong ones are of low confidence.
10. **A second run gives the same output.** The assembly that the folder holds under two names
   was run under both, 43 minutes apart, before the duplicate was noticed. The issue list, the
   copy-number segments, the depth anomalies, the regions that cannot be assessed, the fitted
   model and the scored issues are identical file for file. It is the only repeat with
   unchanged code that was compared.

## Not yet scored

`lc24_nopolish` (ONT only, unpolished) is in `assemblies.tsv` and not in the tables: its run had
not finished when this record was written. `scripts/run_hg002_assemblies.sh` runs whatever is
not done, and `scripts/scoreboard.py` rebuilds the tables and every sum quoted above
(`scoreboard.md`).

## Files

- `assemblies.tsv`: name, FASTA, haplotype patterns and description of each assembly; the
  duplicate is commented out there.
- `details.py` and its three tables for `hprc_year2`: the missed collapses, whole-copy calls by
  length, and the clusters of read sequence the assembly lacks.
- `unconfirmed.py` and `unconfirmed_high.tsv`: every high-confidence copy-number issue the truth
  does not confirm, with the share of it that a truth discordance covers.
- `scoreboard.tsv` (every number per assembly), `scoreboard.by_size.tsv`, `scoreboard.md` (the
  two tables above and the sums quoted in the text, as written by the script),
  `scoreboard.png` with its numbers in `scoreboard.json`.
