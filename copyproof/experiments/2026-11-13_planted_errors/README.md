# Rung 2: errors planted in a truth genome

**Question.** Edit T2T-HG002 Q100 v1.1 so that it carries assembly errors of known class, size and
place, and run copyproof with HG002's real short reads (which come from the unedited genome). Is
each edit found, as the right class, and does anything else change (GUIDE §8 rung 2; decisions D1,
D4, D6)? The first part of this record is from short reads alone. The part
["With long reads"](#with-long-reads-added-2026-11-14) asks what the long-read audit says at
the planted false joins.

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
partners 29 of 34, `MISSING_HAPLOTYPE` 32 of 39; `expected_scored.tsv` lists every expected signal
with the issue that answers it and its confidence). The table above and that file are from the
final classifier (the one that came out of rung 4, see `../2026-11-13_hg002_assemblies/`).

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
   unsupported assembly k-mers, 23 and 29 of them with normally covered flanks
   (`breakpoints.tsv`, written by `scripts/eval_planted.py`). But rung 3 found 1,123 such runs on the unedited genome, nearly all in sequence
   that long reads confirm. So an unsupported run is reported only with a second witness: a
   missing allele from the reads (`BASE_ERROR`), or long reads that also break there
   (`FALSE_JOIN`, needs `--long-read-bam`; the tables above are from short reads alone, so the
   class is silent in them). A join through N has no k-mers at all and is the business of the
   join audit (S4). What the long reads say at these joins was measured a day later: see
   "With long reads" below.
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
     windows number 105 (`../2026-11-13_hg002_q100_false_flags/ambiguous_by_length.tsv`). Such
     segments are now reported, at low confidence.
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

## With long reads (added 2026-11-14)

**Question.** With short reads alone a planted false join is a run of assembly k-mers that no
read holds, and it is not reported. The long-read audit (`--long-read-bam`) is meant to turn such
a run into a `FALSE_JOIN`, and a run of N that no read links into an `UNSUPPORTED_JOIN`. It had
been run only where long reads were aligned to the genome under audit, and none are aligned to
the planted genome. What does it say at the planted false joins when it is given real long reads
(NOTEBOOK entry i, "not done"; the review of 2026-11-14, point 3.4)?

**Reads.** An alignment to the planted genome was made from the public alignment of HiFi reads
to the unedited genome (`hg002_q100_hifi_bam` in `../../resources.tsv`; the one rung 3 reads).

1. The stretches to fetch, in the coordinates of the unedited genome: every edit and the stretch
   homologous to it on the other haplotype, each with 50 kb on both sides. 222 stretches, 33.8 Mb
   (`long_reads.regions.bed`).
2. Every alignment record of the public file in those stretches, fetched by region: 102,425
   records. The 80,615 primary records hold one whole read each: 1.30 Gb of HiFi reads, 16.2 kb
   long on average, 34.6 deep inside the stretches.
3. Those reads, aligned to `plant1.fa` with minimap2 2.30 (`-ax map-hifi`). All of them align
   (`long_reads.reads.tsv`).

The homologous stretches are fetched too because each read is placed again wherever it fits
best: a read of the other haplotype can leave a secondary record at the edit, and a read from
the edit can move to the other haplotype.

**Run.** `scripts/run_planted_hifi.sh DATA`, which does this (`long_read_audit.py` is in this
folder):

    long_read_audit.py regions --edits plant1.edits.tsv --fai hg002v1.1.fasta.gz.fai --flank 50000 --out regions.bed
    samtools view -b -o parts/N.bam URL SEQUENCE:START-END   # one request per stretch, six at a time; samtools cat makes q100.bam
    samtools fastq -F 0x900 q100.bam | gzip > reads.fastq.gz
    minimap2 -ax map-hifi plant1.fa reads.fastq.gz | samtools sort -o plant1.hifi.bam -
    samtools view -b -F 0x100 -o plant1.hifi.no_secondary.bam plant1.hifi.bam
    long_read_audit.py reads --regions regions.bed --fetched q100.bam --aligned plant1.hifi.bam --out run35_hifi.reads.tsv
    scripts/run_sample.py --out run35_hifi --asm mat=plant1.fa --group-map plant1.groups.tsv --db sample=HG002.35x.cpkdb \
        --reads R1 R2 --long-read-bam plant1.hifi.bam
    long_read_audit.py score --edits plant1.edits.tsv --expected plant1.expected.bed --regions regions.bed --run run35_hifi \
        --short-run run35 --truth-run TRUTH_RUN_WITH_LONG_READS --out run35_hifi.audit
    long_read_audit.py records --edits plant1.edits.tsv --expected plant1.expected.bed --bam plant1.hifi.bam --out run35_hifi.audit
    scripts/eval_planted.py --expected plant1.expected.bed --edits plant1.edits.tsv --run run35_hifi --baseline TRUTH_RUN --out run35_hifi.eval

`run35_hifi` holds links to the scan, the windows, the list of unsupported runs and the recruited
reads of `run35`, as the truth genome's run with long reads does, so only the audit is new work.
The last four commands run a second time, as `run35_hifi_no_secondary`, on the alignment without
its secondary records. The fetch took 7 minutes, the alignment 4 (24.5 GB of memory) and each
run of the pipeline about 3. No threshold was changed, before the run or after it.

### Result

A planted false adjacency is the join left by a deletion, either end of an inversion, either run
of N of an inversion that was rejoined through N, and the junction between the two copies of a
tandem duplication. (The last kind was not counted before; it is derived from `edits.tsv`.) The
issue counted is the one within 200 bp of the adjacency, as in `scripts/eval_planted.py`.

| Planted false adjacency | Planted | Shows as a run of unread assembly k-mers | ... which the audit asks the long reads about | `FALSE_JOIN` | `BASE_ERROR` | `UNSUPPORTED_JOIN` | No issue |
|---|---|---|---|---|---|---|---|
| Inversion breakpoint | 24 | 24 | 23 | 23 | 0 | | 1 |
| Join left by a deletion | 36 | 34 | 29 | 21 | 2 | | 13 |
| ... in unique sequence | 20 | 19 | 16 | 13 | 1 | | 6 |
| ... in a segmental duplication | 16 | 15 | 13 | 8 | 1 | | 7 |
| Join through 100 N (two per inversion) | 8 | 0 | (join audit: 8) | | | 8 | 0 |
| Junction of a tandem duplication | 39 | 38 | 7 | 5 | 2 | | 32 |

All 49 `FALSE_JOIN` issues are of high confidence; the 8 `UNSUPPORTED_JOIN` are of low
confidence (`long_reads.summary.tsv`; one row per adjacency in `long_reads.adjacencies.tsv`).

`FALSE_JOIN` by size of the edit, of all planted and, in brackets, of those the audit asked about
(`long_reads.by_size.tsv`):

| | 2 kb | 5 kb | 10 kb | 50 kb | 200 kb |
|---|---|---|---|---|---|
| Inversion breakpoint | | | 7/8 (7/7) | 8/8 (8/8) | 8/8 (8/8) |
| Join left by a deletion | 0/8 (0/5) | 5/8 (5/8) | 7/8 (7/7) | 5/8 (5/5) | 4/4 (4/4) |
| Junction of a tandem duplication | 0/8 (0/0) | 0/8 (0/1) | 0/8 (0/1) | 3/8 (3/3) | 2/7 (2/2) |

**Nothing is reported away from the edits.** The audit made or changed 86 issues: the 49
`FALSE_JOIN`, the 8 `UNSUPPORTED_JOIN`, 4 `BASE_ERROR` and 25 `EXPANSION` (planted duplications
whose unread k-mers long reads do not hold either; 7 of them rose from low to high confidence).
Every one lies within 5 kb of an edit (`long_reads.audit_issues.tsv`, `long_reads.issue_changes.tsv`).
The issue list grows from 656 rows to 717; `scripts/eval_planted.py` counts 298 at planted edits
where it counted 237, and 419 elsewhere as before. Its table of rung 2, written again from this
run (`long_reads.sensitivity.tsv`), differs from `sensitivity.tsv` in the nine rows for the join
of a deletion and in no other: 23 of 36 joins have a `FALSE_JOIN` or `BASE_ERROR`, where none had.

Inside the fetched stretches the audit asked about 4 runs that lie more than 5 kb from every
edit. Long reads hold all their k-mers, as they do in the truth genome's run with the public
alignment, and none becomes an issue (`long_reads.elsewhere.tsv`). Four is a small number: see
"What this test cannot show".

**The same audit without the alignment's secondary records** (`long_reads_no_secondary.*`):

| | Alignment as minimap2 writes it | Secondary records removed |
|---|---|---|
| Secondary records | 85,778 | 0 |
| Inversion breakpoints reported as `FALSE_JOIN` | 23 of 24 (23 high) | 8 of 24 (5 high) |
| Joins left by a deletion reported as `FALSE_JOIN` | 21 of 36 (21 high) | 3 of 36 (2 high) |
| Junctions of a tandem duplication reported as `FALSE_JOIN` | 5 of 39 (5 high) | 0 of 39 |
| Joins through N reported as `UNSUPPORTED_JOIN` | 8 of 8 (low) | 8 of 8 (low) |
| Records clipped at an inversion breakpoint, median of those asked about | 27 | 2 |
| Records clipped at the join of a deletion, median | 23 | 1 |

### What was learned

1. **With a minimap2 alignment that keeps its secondary records, the audit reports planted false
   joins as `FALSE_JOIN`, and reports nothing away from the edits.** 23 of 24 inversion
   breakpoints; of the joins left by deletions, 16 of 20 from 10 kb, 5 of 8 at 5 kb and none of 8
   at 2 kb. All are of high confidence. Each of the 49 `FALSE_JOIN` issues lies at a planted
   adjacency (5 are junctions of tandem duplications), and every issue the audit made or changed
   lies at an edit.
2. **A join is missed when the long reads are not asked, or when alignments run through it.**
   - *Not asked.* The audit asks about a run of unread assembly k-mers only if short reads cover
     its flanks at 0.85 copies or more. That leaves out 1 inversion breakpoint, 5 deletion joins
     and 31 duplication junctions; 2 more deletion joins and 1 duplication junction show no run
     of unread k-mers at all. Every deletion join of 10 kb or more that the audit asked about is
     a `FALSE_JOIN` (16 of 16).
   - *Run through.* At the 8 joins of 2-kb deletions 405 alignments run through the join, and
     every one of them holds an insertion of about 2,000 bases there: the deleted sequence. 402
     are secondary records of reads whose primary record is on the other haplotype, and none has
     a mapping quality of 20 or more (`long_reads.records_by_size.tsv`). The audit takes an
     alignment through a run as a read through it, so none of the 5 joins it asked about is a
     `FALSE_JOIN`. Two are reported as `BASE_ERROR` ("a wrong base or small consensus error"),
     which names a missing 2 kb badly, and three are not reported. At 5 kb, 6 such alignments
     remain in all, and they cost three of the eight joins: one has 3 alignments through, and two
     have 1 or 2, which the rule counts neither as a read through nor as none.
3. **The call comes from alignment, not from the reads' content, and from secondary records.**
   The audit's first question is whether the long reads over a run hold its unread k-mers. At
   the 59 planted adjacencies it asked about, no long read holds any of them. But the verdict
   "not carried" needs 5 reads with their sequence over the run, and the median is 2, in
   stretches where the reads lie 34.6 deep: a read that crosses a false join fits the other
   haplotype better, its primary record goes there, and what stays at the join is a secondary
   record, which holds no sequence. So content says "not carried" at 8 of the 59 and "unclear"
   at 51, and the second rule decides: no alignment through the run and at least 3 clipped at
   it. Of the 582 records clipped at the 24 inversion breakpoints, 534 are secondary, and all 534
   belong to reads whose primary record is on the other haplotype; 34 are primary. At the 36
   deletion joins the same counts are 857, 815, 815 and 41. Take the secondary records away and
   the median number of clipped records at a join the audit asks about falls from 27 to 2
   (inversions) and from 23 to 1 (deletions), and `FALSE_JOIN` is called at 8 of 24 and 3 of 36.
   An alignment file without secondary records loses most of this class.
4. **A join through N is reported, but with the confidence of an ordinary gap.** No read links
   either side of any of the 8 runs of N (27 to 45 reads on each side, 8 to 17 of them clipped
   at the run), so all 8 are `UNSUPPORTED_JOIN`. They are of low confidence because an unlinked
   join is made firm only by another issue within 20 kb (in a real release 34 of 42 scaffold
   joins have no linking read, rung 1), and an inversion leaves no other issue beside it. So an
   inversion that was closed with gaps cannot be told from a scaffold gap by this channel.
5. **The junction of a tandem duplication is seldom asked about, and the duplication is found
   anyway.** The copies on both sides of the junction are short of reads themselves, so only 7 of
   39 junctions pass the flank rule: 5 are `FALSE_JOIN` (50 and 200 kb) and 2 are `BASE_ERROR`
   (5 and 10 kb). At duplications of 2 to 10 kb, 905 alignments run through the junction, 397 of
   them primary, each with a deletion of about the copy's length. The copy-number channel finds these
   duplications from 5 kb (first table of this record), and the long reads raise 7 of the
   `EXPANSION` issues from low to high confidence.
6. **Two aligners place the same reads differently.** minimap2 2.30 put the primary record of
   9,342 of the 80,615 reads on another sequence than Winnowmap 2.03 had, gave 45,373 of them a
   mapping quality of 20 or more where Winnowmap gave 65,119, and wrote 85,778 secondary records
   where Winnowmap wrote 21,697 (`long_reads.reads.tsv`). Since the calls rest on secondary
   records, the numbers above are those of this aligner.

### What this test cannot show

- **The reads were chosen by their alignment to the unedited genome.** A read was aligned again
  only if its primary record in the public file lies in a fetched stretch. Reads that an
  alignment of the whole read set to the planted genome would bring to an edit from elsewhere
  (a paralogous copy, for the edits in segmental duplications) are missing. 645 reads have a
  secondary or supplementary record in the stretches and their primary record outside them;
  reads with no record there cannot be counted.
- **Only the neighbourhood of the edits is covered:** 33.8 Mb of 6 Gb. Everywhere else the audit
  had no reads and said nothing, so "nothing away from the edits" is a statement about the 4
  runs it asked about in the fetched stretches. The count of `FALSE_JOIN` on a whole correct
  genome is the one of rung 3, with the public alignment: 14, of which 13 on known issues.
- **The planted joins are clean, and the other haplotype is right at each of them.** Every
  inversion is in unique sequence. A false join inside a repeat, one that both haplotypes share,
  or one in an assembly with a single haplotype was not planted. There the reads would have
  nowhere better to go, and the audit would see other records than it saw here.
- **One read set and one aligner:** HiFi reads at a depth of 34.6 per haplotype, minimap2 2.30 with
  `-ax map-hifi`. The public alignments of this genome were made with Winnowmap, which wrote a
  quarter as many secondary records for the same reads in the same stretches (21,697). What the
  audit says on a Winnowmap alignment of the planted genome, or with ONT reads, was not run.
- **Whether these numbers hold with other thresholds.** None was changed or varied here.

## Files

- `sensitivity.tsv`: found and planted per kind, signal, stratum and size.
- `expected_scored.tsv`: the 266 expected signals one by one: found or not, the class and role of
  the issue that answers it, and that issue's confidence (written by `scripts/eval_planted.py`,
  as is `sensitivity.tsv`, which counts its rows).
- `breakpoints.tsv`: per planted false adjacency, whether it shows as a run of unsupported
  assembly k-mers, and whether that run has normally covered flanks.
- `edits.tsv`: every planted edit in original and edited coordinates.
- `new_elsewhere.tsv`: the 13 issues that are neither at an edit, nor at its homologue, nor in the
  unedited genome.

With long reads (all written by `long_read_audit.py` in this folder, except the last):

- `long_reads.regions.bed`: the 222 stretches of the unedited genome whose reads were fetched
  (`regions`).
- `long_reads.reads.tsv`: what was fetched and what the new alignment made of it (`reads`).
- `long_reads.adjacencies.tsv`: one row per planted false adjacency: the run of unread k-mers,
  whether the audit asked about it, what it counted, both verdicts and the issue (`score`).
- `long_reads.summary.tsv`, `long_reads.by_size.tsv`: the same counted per kind of adjacency and
  stratum, and per kind and size of the edit.
- `long_reads.audit_issues.tsv`: every issue the audit made or changed, with its distance to
  the nearest edit.
- `long_reads.elsewhere.tsv`: the runs the audit asked about inside the fetched stretches and
  away from the edits, beside the verdicts of the truth genome's run at the same places.
- `long_reads.issue_changes.tsv`: issues by class and confidence, without and with long reads.
- `long_reads.records.tsv`, `long_reads.records_by_size.tsv`: the alignment records at each
  planted adjacency by type: through it, clipped at it, with their sequence over it (`records`).
- `long_reads_no_secondary.*`: the same tables for the run without secondary records.
- `long_reads.sensitivity.tsv`: `sensitivity.tsv` again, from the run with long reads
  (`scripts/eval_planted.py`).

To write them again from the run directories (`P` is `DATA/work/plant`, `X` this folder):

    long_read_audit.py regions --edits $P/plant1.edits.tsv --fai DATA/hg002/asm/hg002v1.1.fasta.gz.fai --out $X/long_reads.regions.bed
    long_read_audit.py reads --regions $X/long_reads.regions.bed --fetched $P/hifi/q100.bam --aligned $P/hifi/plant1.hifi.bam --out $X/long_reads.reads.tsv
    long_read_audit.py score --edits $P/plant1.edits.tsv --expected $P/plant1.expected.bed --regions $X/long_reads.regions.bed \
        --run $P/run35_hifi --short-run $P/run35 --truth-run DATA/work/hg002_v1.1/run35_hifi --out $X/long_reads
    long_read_audit.py records --edits $P/plant1.edits.tsv --expected $P/plant1.expected.bed --bam $P/hifi/plant1.hifi.bam --out $X/long_reads
    cp $P/run35_hifi.eval.sensitivity.tsv $X/long_reads.sensitivity.tsv

and `score` and `records` again with `run35_hifi_no_secondary`, `plant1.hifi.no_secondary.bam`
and `--out $X/long_reads_no_secondary`.
