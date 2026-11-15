# What the missing alleles in neither parent's reads are: the release's own long reads, and the changed bases

**Question.** The record on missing alleles by the parents' reads
(`../2026-11-14_missing_alleles_by_parent/README.md`) found that 545 of 1,701 and 922 of 1,773
strong missing-allele placements in the HPRC releases of HG01943 and NA19240 are in neither
parent's reads, and could not say what they are. Three readings were open:

1. the release's own reads lack the allele: the DNA of the short reads and the DNA of the release
   differ;
2. the release's own reads hold both alleles: a change in part of the cells;
3. the release's own reads hold the allele and the release does not: an error of the assembly.

Two sources that the earlier record did not use are read here: the long reads each release was
built from, and the changed bases themselves.

No rule was made and no code under `src/` or `engine/` changed. Every sample here had been looked
at before, so nothing is held out.

**Inputs.** Four runs of the earlier record, each with the child's short reads and both parents'
(`runs.tsv`):

| Run | Strong missing-allele placements | In a parent's reads | In neither | Events: in a parent, in neither | Long reads |
|---|---|---|---|---|---|
| NA19240, HPRC release 2 | 1,773 | 851 | 922 | 467, 460 | fetched |
| HG01943, HPRC release 2 | 1,701 | 1,156 | 545 | 731, 276 | fetched |
| HG002, `hprc_year2` | 996 | 925 | 71 | 469, 38 | not looked at |
| HG002, truth genome | 74 | 40 | 34 | 26, 18 | not looked at |

The placements in a parent's reads are the comparison throughout: inherited alleles that the
release lacks.

The long reads are HPRC's own: for each release it publishes the HiFi reads the release was built
from, aligned to that release and indexed. The addresses are the `BAM` variables of
`scripts/run_na19240_heldout.sh` and `scripts/run_hg01943.sh`, and nothing was fetched from
anywhere else. Both answered on 14 November 2026. By their headers the files were made with
minimap2 2.28 (`-x lr:hqae`) from 3 read files for NA19240 and 5 for HG01943, and they list every
sequence of both haplotypes, 129 and 244: the reads were aligned to the two haplotypes together.
Fetched were the records at the placements only, kept without base qualities:

| | Stretches asked for | Records | Read sequence | Kept as | Time |
|---|---|---|---|---|---|
| NA19240 | 1,744 | 23,079 | 434 Mb | 104 MB | 15 min |
| HG01943 | 1,603 | 39,720 | 777 Mb | 189 MB | 16 min |

**Run.**

    .venv/bin/python experiments/2026-11-14_missing_alleles_long_reads/fetch_long_reads.py DATA
    .venv/bin/python experiments/2026-11-14_missing_alleles_long_reads/tables.py DATA experiments/2026-11-14_missing_alleles_long_reads

`DATA` is the data root. The first command needs the network and writes the reads to
`DATA/work/missing_alleles_long_reads`. The second needs none, takes 3 minutes and writes every
table here. It also writes one row per
placement to `DATA/work/missing_alleles_long_reads/placements.tsv` (1.4 MB, so not in this folder).

**Terms.** As in the earlier record: a *placement* is one row of a run's `missing_clusters.tsv`,
an *event* is the placements that hold mostly the same reads (`twin_group`; usually the two
haplotypes' places of one allele), *strong* is 6 reads or more and 5 tracked k-mers or more, *in a
parent's reads* means half or more of the placement's tracked k-mers are in the mother's or in the
father's read database, and the *read count* is the median count of those k-mers in the child's
short reads. New here: a long read *holds the allele* or *holds the release's sequence* at a
placement as defined in the next section.

## How a placement is sorted from its long reads

This was written into NOTEBOOK entry u at 15:45, before any long read was searched for a k-mer, and
`placements.py` holds it as code. By then the placements' k-mers had been recovered and one
position had been fetched to measure time and size. The times are those of one machine and
nothing was pushed, so the order rests on the notebook and on `fixed_beforehand.sha256`.

- Long reads: the records of the two alignments that are neither secondary nor supplementary and
  have their sequence.
- A read *covers* a placement when its alignment covers the placement's interval, from the base
  before `start` to `end`, and 100 bases on either side. Only the stretch of the read aligned to
  that window is searched.
- The *allele* is the placement's tracked k-mers. The *release's sequence there* is the tracked
  k-mers of the release that overlap the interval and that none of the placement's short reads
  holds. A read holds one or the other when its stretch holds one or more of those k-mers, on
  either strand.
- The *share* of a placement is the reads that hold the allele over the reads that hold the allele
  or the release's sequence.
- Groups: *too few reads to say* (fewer than 10 reads that hold either); *no long read holds the
  allele*; *a minority do* (share under 0.25); *about half do* (0.25 to 0.75); *nearly all do*
  (above 0.75). With 10 reads, none with the allele is seen with a chance under 5% when the true
  share is 0.26 or more.
- An *event* is sorted from every row of the cluster table with its `twin_group`: it takes the
  highest group among them in the order nearly all, about half, a minority; it counts as "no long
  read holds the allele" only when each of its rows has 10 reads or more and none holds the allele.
- Readings: "the release's own reads lack the allele" is the group with no long read; "they hold
  both alleles" is a minority or about half, with reads that hold the release's sequence beside
  them; "they hold the allele and the release does not" is nearly all.
- Known beforehand to blur this: where the release holds a stretch once, or its two haplotypes are
  alike over the length of a read, reads of both haplotypes lie at one place and an inherited
  allele gives about half there.

Two things the block did not settle and the script had to: the window ends 100 bases after the base
at `end`, counting that base; and a placement where the release has no tracked k-mer of the kind
described (10 of 1,773 on NA19240, 8 of 1,701 on HG01943) is counted as too few to say.

The tables that follow from the block are `long_reads_events.tsv`, `long_reads_by_parent.tsv`,
`long_reads_by_read_count.tsv`, `long_reads_by_other_assemblies.tsv`,
`long_reads_share_steps.tsv` and `long_reads_aligned_base.tsv`. The others were added after the
counts had been seen, among them the two that the reading of "both alleles" leans on:
`long_reads_event_patterns.tsv` and `share_of_all_reads.tsv`.

## Result

**The k-mers of every placement were recovered** (`checks.tsv`). A cluster table does not store
which reads and k-mers made a placement. They were collected again from each run's `missing/`
folder with the package's own `classify_reads` and the loop of `cluster_reads`. All 10,953 rows of
the four cluster tables come back under their cluster number with the same sequence, kind, reads,
strands, tracked k-mers, read count, start and parental shares; so do all 4,544 strong
missing-allele rows.

**The long reads are deep enough at placements in neither parent's reads, and often not at those
in a parent's on NA19240** (`long_reads_reads.tsv`). A placement in neither parent's reads is
covered by a median of 17 long reads on NA19240 and 32 on HG01943. A placement in a parent's reads
is covered by a median of 8 and 33: on NA19240, 469 of those 851 placements have fewer than 10
covering reads. Of the 33,927 covering reads at placements in neither parent's reads, 437 hold
neither the allele nor the release's sequence and 32 hold both.

**At nearly every event in neither parent's reads, a part of the long reads at one place hold the
allele** (`long_reads_events.tsv`). Events, sorted from all their rows:

| | Events | No long read holds the allele | A minority do | About half do | Nearly all do | Too few reads |
|---|---|---|---|---|---|---|
| NA19240, in neither parent's reads | 460 | 2 | 42 | 399 | 2 | 15 |
| HG01943, in neither parent's reads | 276 | 1 | 31 | 244 | 0 | 0 |
| NA19240, in a parent's reads | 467 | 24 | 19 | 68 | 79 | 277 |
| HG01943, in a parent's reads | 731 | 15 | 3 | 458 | 227 | 28 |

Placements, each by itself (`long_reads_by_parent.tsv`):

| | Placements | No long read holds the allele | A minority do | About half do | Nearly all do | Too few reads |
|---|---|---|---|---|---|---|
| NA19240, in neither parent's reads | 922 | 420 | 63 | 400 | 2 | 37 |
| HG01943, in neither parent's reads | 545 | 237 | 64 | 244 | 0 | 0 |
| NA19240, in a parent's reads | 851 | 141 | 59 | 98 | 79 | 474 |
| HG01943, in a parent's reads | 1,156 | 159 | 130 | 582 | 235 | 50 |

The placements with no long read are mostly the other haplotype's place of an event that has its
allele at the first. Where about half hold the allele, the others hold the release's sequence:
summed over those placements in neither parent's reads, 3,077 of 7,086 long reads hold the allele
on NA19240 (43%) and 3,120 of 7,967 on HG01943 (39%). At the placements with no long read, 7,323
and 7,705 reads hold the release's sequence and none the allele.

**The two places of one event differ, and not as they do for inherited alleles**
(`long_reads_event_patterns.tsv`; made after the counts were seen). Most events have two rows, one
on each haplotype. How those two fall together:

| Events with two rows | NA19240, in neither | HG01943, in neither | NA19240, in a parent | HG01943, in a parent |
|---|---|---|---|---|
| about half at one, none at the other | 357 | 214 | 11 | 11 |
| a minority at one, none at the other | 32 | 16 | 1 | 1 |
| nearly all at one, none or a minority at the other | 2 | 0 | 73 | 199 |
| about half at both | 1 | 0 | 29 | 110 |
| none at both | 2 | 1 | 19 | 12 |
| all events with two rows | 447 | 262 | 373 | 416 |

An inherited allele that the release lacks shows in one of two ways. Where the long reads of the
two haplotypes are told apart, nearly all at one place hold it (73 and 199 events). Where they are
not, about half hold it at both places (29 and 110), or at the only place when the event has one
row (13 and 290 events): that is the blur named beforehand. The alleles in neither parent's reads
show a third pattern, which 11 and 11 inherited events show: the haplotypes are told apart, since
the other place holds no read with the allele, and still only a part of the reads at the allele's
own place hold it. The two kinds of place are equally deep (`long_reads_by_parent.tsv`: a median
of 17 reads that hold either at the placements with no long read and 17 at those where about half
hold the allele on NA19240, 32 and 32 on HG01943). A second copy of the stretch at the allele's
place would bring its own reads, so it does not explain the pattern.

**Of all the long reads of an event, a fifth hold the allele; of an inherited allele's, a third to
a half** (`share_of_all_reads.tsv`; made after the counts were seen). Summed over all rows of
an event, for events with 20 long reads or more:

| | Events | Share of the long reads with the allele: summed, median, quartiles | Read count of the short reads over two copies: median, quartiles |
|---|---|---|---|
| NA19240, in neither parent's reads | 452 | 0.206 (3,271 of 15,901); 0.206; 0.151 to 0.257 | 0.293; 0.275 to 0.367 |
| HG01943, in neither parent's reads | 276 | 0.189 (3,387 of 17,927); 0.190; 0.156 to 0.226 | 0.281; 0.259 to 0.316 |
| NA19240, in a parent's reads | 208 | 0.321 (2,341 of 7,304); 0.373; 0.250 to 0.444 | 0.504; 0.403 to 0.587 |
| HG01943, in a parent's reads | 699 | 0.457 (19,644 of 42,990); 0.467; 0.418 to 0.509 | 0.509; 0.422 to 0.597 |

The short reads' figure cannot go below a quarter: a k-mer is placed from 7 reads (the earlier
record's `runs.tsv`), and two copies give 27.3 and 28.5.

**Neither the read count nor the other assemblies sort these events**
(`long_reads_events_by.tsv`, `long_reads_by_read_count.tsv`, `long_reads_by_other_assemblies.tsv`).
Events in neither parent's reads where a minority or about half of the long reads hold the allele:

| | NA19240 | HG01943 |
|---|---|---|
| read count under 9 | 234 of 245 | 162 of 162 |
| read count 9 to 12.5 | 195 of 203 | 108 of 109 |
| read count 13 or more | 12 of 12 | 5 of 5 |
| both HGSVC3 assemblies hold the allele | 173 of 180 | |
| the Verkko assembly only | 67 of 68 | |
| the hifiasm assembly only | 61 of 63 | |
| neither HGSVC3 assembly holds it | 140 of 149 | |

**The second count agrees** (`long_reads_aligned_base.tsv`, `long_reads_by_aligned_base.tsv`).
Where the change is one base, the base each covering read has at that position by its alignment
was counted beside its k-mers. At placements in neither parent's reads, 6,182 of 6,213 reads that
hold k-mers of the allele have the reads' base there, and 24,311 of 24,314 reads that hold only
k-mers of the release have the release's base. Sorted by the aligned base in place of the k-mers,
816 of 826 placements on NA19240 and 500 of 506 on HG01943 fall in the same group.

**Nearly all placements in neither parent's reads are one changed base** (`changed_bases.tsv`,
`changed_bases_detail.tsv`). The release's sequence at each placement was taken with
`cpk extract --regions`, and each tracked k-mer of the placement was laid on it without gaps where
it fits best. A placement is one changed base when every one of its k-mers differs from the
release at one and the same position by one and the same base.

| | In neither parent's reads: placements | events | In a parent's reads: placements | events |
|---|---|---|---|---|
| NA19240 | 826 of 922 | 425 of 460 | 354 of 851 | 207 of 467 |
| HG01943 | 506 of 545 | 261 of 276 | 583 of 1,156 | 409 of 731 |
| HG002, `hprc_year2` | 70 of 71 | 38 of 38 | 399 of 925 | 207 of 469 |
| HG002, truth genome | 32 of 34 | 17 of 18 | 23 of 40 | 15 of 26 |

An event is counted once, under the smallest difference among its placements. An allele belongs
to one haplotype; at the other haplotype's place it differs from the release by its own change and
by whatever differs between the two haplotypes there. The long reads agree with that
(`long_reads_at_the_nearer_placement.tsv`): of the events in neither parent's reads with one
placement of one changed base and one of a larger difference, the long reads hold the allele at
the placement of one changed base alone in 16 of 20 on NA19240 and 7 of 7 on HG01943, and at the
other alone in 1.

The placements that could not be resolved to one changed base, and why:

| Placements | In neither: NA19240, HG01943, `hprc_year2`, truth genome | In a parent's: the same four |
|---|---|---|
| two or more changed bases less than 31 apart | 44, 16, 1, 0 | 58, 62, 48, 5 |
| the reads have 1 base more or fewer | 25, 14, 0, 0 | 186, 295, 264, 0 |
| the reads have 2 to 30 bases more or fewer | 14, 4, 0, 2 | 73, 64, 153, 0 |
| the reads have more than 30 bases more or fewer | 0, 0, 0, 0 | 2, 4, 6, 0 |
| not determined | 13, 5, 0, 0 | 178, 148, 55, 12 |

"Not determined" is: a k-mer fits the release at more than one place (44 placements in all), the
k-mers hold two different bases at one position (7), or a k-mer does not fit without gaps and the
placement's reads, laid over one another, give no single difference between two anchors in the
release (360). Counted by events, two or more changed bases are 12, 4, 0 and 0 events in neither
parent's reads.

A second way of deriving the change was run beside the first (`checks.tsv`): the placement's reads
laid over one another and compared with the release between two anchors. Of the 2,793 placements of
one changed base, it gives the same base at the same position at 2,686 (at 60 of them beside other
differences further out), another length between the anchors at 21, another base at none, and no
anchored sequence at 86.

**The cluster table's `size` of 0 does not mean the same length** (`size_column.tsv`). It is 0
whenever one chain of anchors runs through the missing k-mers, and a chain joins anchors whose
diagonals lie up to 30 bases apart (`diag_tol` in `src/copyproof/missing.py`), so a length change
of up to 30 bases stays inside one chain. No resolved placement with `size` 0 has a longer one. Of the
placements with `size` 0, the reads have another length at 39 of 919 on NA19240,
18 of 543 on HG01943, 0 of 71 on `hprc_year2` and 2 of 34 on the truth genome when they are in
neither parent's reads, and at 256 of 814, 359 of 1,128, 415 of 911 and 0 of 36 when they are in a
parent's.

**The changed bases are of other kinds than the inherited ones** (`substitutions.tsv`). Events of
one changed base, written from the pyrimidine, from the release's base to the reads':

| | Events | C>A | C>G | C>T | T>A | T>C | T>G | C>T or T>C |
|---|---|---|---|---|---|---|---|---|
| NA19240, in neither parent's reads | 425 | 70 | 20 | 104 | 42 | 98 | 91 | 202 (48%) |
| HG01943, in neither parent's reads | 261 | 27 | 12 | 85 | 19 | 59 | 59 | 144 (55%) |
| HG002 `hprc_year2`, in neither | 38 | 9 | 4 | 12 | 4 | 7 | 2 | 19 (50%) |
| HG002 truth genome, in neither | 17 | 4 | 1 | 6 | 2 | 3 | 1 | 9 (53%) |
| NA19240, in a parent's reads | 200 | 14 | 28 | 59 | 22 | 63 | 14 | 122 (61%) |
| HG01943, in a parent's reads | 403 | 41 | 44 | 120 | 21 | 139 | 38 | 259 (64%) |
| HG002 `hprc_year2`, in a parent's | 206 | 15 | 15 | 80 | 19 | 57 | 20 | 137 (67%) |
| HG002 truth genome, in a parent's | 15 | 0 | 2 | 3 | 1 | 5 | 4 | 8 (53%) |

In the two releases T>G is 21% and 23% of the events in neither parent's reads (91 of 425, 59 of
261) and 7% and 9% of those in a parent's (14 of 200, 38 of 403). With the neighbouring bases
(`substitutions_with_neighbours.tsv`): the T>G events in neither parent's reads are followed by A
or T in 75 of 91 and 52 of 59, those in a parent's reads in 7 of 14 and 17 of 38. C>T stands before
a G in 41 of 104 and 11 of 85 events in neither parent's reads, and in 19 of 59 and 38 of 120 in a
parent's. For an inherited allele the direction is not known: the release may hold the newer base.

The kinds do not differ in what else the tables say of them
(`substitutions_beside_the_cluster_table.tsv`, `substitutions_by_long_reads.tsv`). In the two
releases every kind in neither parent's reads has a median read count of 8 to 9 and 35% to 41% of
its reads on the rarer strand, and the long reads hold a T>G allele about as often as any other:
at 86 of 91 events on NA19240 and 58 of 59 on HG01943, a minority or about half do.

## The three readings

Events in neither parent's reads, from `long_reads_events.tsv`:

| Reading | Long reads | NA19240 | HG01943 |
|---|---|---|---|
| the release's own reads lack the allele | none holds it, 10 reads or more at each place | 2 of 460 | 1 of 276 |
| they hold both alleles | a minority or about half hold it, the others the release's sequence | 441 of 460 | 275 of 276 |
| they hold the allele and the release does not | nearly all hold it | 2 of 460 | 0 of 276 |
| not said | too few reads | 15 of 460 | 0 of 276 |

The tables separate the three. The first and the third are 5 of 736 events. The second is 716.

Whether "both alleles" is a change in part of the cells was not sorted by a rule fixed beforehand,
and rests on the two tables made afterwards: the allele's reads lie at one haplotype's place and
are a part of the reads there (about half at one place and none at the other in 571 of the 709
events with two rows, a minority and none in 48 more), where inherited alleles are nearly all of
one place or about half of both; the two places are equally deep; and the allele is in a fifth of
all long reads of an event, where an inherited allele is in a third to a half. A haplotype
whose molecules come in two kinds, in DNA of one person, is what a change in part of the cells
gives. Cells were not counted.

## What was learned

1. **The release's own reads hold these alleles.** At 716 of the 721 events in neither parent's
   reads that have enough long reads, a part of the reads the release was built from hold the
   allele. The DNA of the short reads and the DNA of the release do not differ in them.
2. **They hold them beside the release's sequence, on one haplotype, in two reads of five.** That
   is the picture of a change in part of the cells, and not the picture inherited alleles give in
   the same alignments.
3. **By its own reads the release is not in error at these places.** It holds the allele that the
   larger part of its reads at that place hold (4,009 of 7,086 and 4,847 of 7,967 where about half
   hold the other). "The reads carry an allele that no sequence in the release holds" stays true of
   such a placement; "the assembler erred" does not follow. Nearly all long reads hold the allele
   at 2 of 736 events.
4. **For inherited alleles the same count does point at the release.** Where the long reads are
   told apart by haplotype, nearly all at the place hold the allele and not the release's sequence
   (79 events on NA19240, 227 on HG01943).
5. **These alleles are single changed bases** (92%, 95%, 100% and 94% of the events), of other
   kinds than the inherited ones, with T>G before A or T standing out in both releases.
6. **What another assembly holds says little about them.** Whether both, one or neither HGSVC3
   assembly holds the allele, the release's long reads hold it in part at 94% to 99% of the events.
   An assembler meets two alleles there and keeps one.
7. **`size` is not a length.** The earlier record read `size` 0 as "the same length as the release".
   Among placements with `size` 0 the reads have another length at 39 of 919 and 18 of 543 in neither
   parent's reads, and at 256 of 814, 359 of 1,128 and 415 of 911 in a parent's reads, where a single
   base more or fewer is the most frequent other case.

## What this does not show

- HG002. No long reads were looked at for its 38 and 18 events; only their changed bases were
  derived.
- Where the change arose. Both kinds of reads come from DNA of a cell line. Whether the changed
  cells were already in the person or arose in culture is not in these tables.
- What made the changes. A signature of T>G changes is described for B-cell cancers whose
  immunoglobulin genes are hypermutated, and is put down, tentatively and in part, to polymerase
  eta in the somatic hypermutation of lymphoid cells (signature 9 of Alexandrov et al., Nature
  500:415-421, 2013, https://doi.org/10.1038/nature12477; SBS9 in COSMIC,
  https://cancer.sanger.ac.uk/signatures/sbs/sbs9/). These cell lines are made from B cells.
  Whether the changes here are of that signature was not tested. The paper names A or T before
  the T; here A or T follows the T in 75 of 91 and 52 of 59 events and stands before it in 57 and
  28 (`substitutions_with_neighbours.tsv`), and the changes were not compared with the
  signature's 96 classes.
- How many such changes the release itself holds. A placement exists only where the release kept
  the unchanged allele. Where it kept the changed one, the release holds sequence that neither
  parent has, and this channel does not see it.
- The true shares. What is counted is cut at both ends: from below by the 7 short reads a k-mer
  needs, from above by the release keeping the other allele.
- Shares under about a quarter at single placements of NA19240, where a place has a median of 17
  reads. Three events in all have no long read with the allele, and 15 have too few reads.
- How few long reads some events rest on (added after the record was checked a second time).
  Where 10 reads or more hold the allele or the release's sequence, one read with the allele sorts
  a placement as "a minority do". Of NA19240's 42 events in neither parent's reads that are sorted
  so, 13 have one or two long reads with the allele at the place that has most of them (7 have
  one, 6 two; counting all places of an event together, 12 have one or two). So 13 of NA19240's
  441 events stand on one or two reads. Each of HG01943's 31 has three or more. These counts are
  from the row-by-row tables that `tables.py` writes to the data root (`placements.tsv` and
  `other_rows_of_the_events.tsv`); no table in this folder holds them.
- Much about inherited alleles on NA19240: 474 of her 851 placements in a parent's reads have too
  few long reads.
- The weak placements, of which 197 on HG01943 and 298 on NA19240 are in neither parent's reads
  (the earlier record).
- A rule. None was made. HG002, HG01943 and NA19240 have all been looked at for this question, so
  a rule would need another person to test it; and a rule that uses long reads needs the long
  reads.
- An order of events that an outside reader can check. The sorting was fixed before the counts by
  this machine's clock only.

## Elsewhere in the repository

- `../../README.md` said of these calls that what they are was not measured. It now gives the
  counts above, and `../README.md` has the row for them.
- The earlier record said, and NOTEBOOK entry t says, that the reads have the same length as the
  release at 543 of 545 and 919 of 922 placements in neither parent's reads. Those are counts of
  `size` 0. By the reads, 522 and 869 of them have the same length, 18 and 39 another, and 3 and
  11 are not determined (`size_column.tsv`). The earlier record's bullet and its table column were
  corrected after this record was written, and entry t carries a note.
- `../../docs/method.md` said that a `MISSING_ALLELE` comes with the difference in length. The
  cluster table gives one only from 31 bases, and that sentence was corrected after this record
  was written.

## Files

- `fetch_long_reads.py`: fetches the long reads at the placements from the two public alignments.
- `placements.py`: recovers the reads and k-mers of every placement, derives the change, counts the
  long reads and holds the sorting rule.
- `tables.py`: writes every table below from the run directories and the fetched reads.
- `fixed_beforehand.sha256`: checksum and time of the sorting block of NOTEBOOK entry u.
- `runs.tsv`, `checks.tsv`: the runs, what was fetched, and the checks of the recovery and of the
  changed bases.
- `long_reads_by_parent.tsv`, `long_reads_events.tsv`, `long_reads_event_patterns.tsv`: placements
  and events per long-read group; which groups the rows of one event fall into together.
- `long_reads_by_read_count.tsv`, `long_reads_by_other_assemblies.tsv`, `long_reads_events_by.tsv`:
  the groups by read count and by what the HGSVC3 assemblies hold, by placements and by events.
- `long_reads_share_steps.tsv`, `share_of_all_reads.tsv`: the share in finer steps; the allele's
  share of all long reads of an event beside the short reads' count.
- `long_reads_reads.tsv`, `long_reads_aligned_base.tsv`, `long_reads_by_aligned_base.tsv`: the
  reads themselves, and the second count by aligned base.
- `changed_bases.tsv`, `changed_bases_detail.tsv`, `size_column.tsv`: what the placements are, with
  every reason; what `size` 0 covers.
- `long_reads_by_change.tsv`, `long_reads_at_the_nearer_placement.tsv`: the long reads by the kind
  of change.
- `substitutions.tsv`, `substitutions_with_neighbours.tsv`, `substitutions_by_long_reads.tsv`,
  `substitutions_beside_the_cluster_table.tsv`: the six kinds of substitution, alone, with the base
  before and after, by long-read group, and beside read count, strands and other assemblies.
