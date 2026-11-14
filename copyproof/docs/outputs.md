# copyproof outputs

`copyproof run --out DIR ...` writes the files below. Coordinates are 0-based and half-open, on the
sequences of the release itself (not on a reference). A flagged locus is an *issue*: a hypothesis
with evidence, to be adjudicated.

## Issues

**`issues.tsv`**: one row per issue.

| Column | Meaning |
|---|---|
| `id` | `cp000001`, ... in order of position |
| `contig`, `start`, `end` | the locus |
| `class` | see the table of classes below |
| `role` | `locus`: the evidence points here. `partner`: copy number is off here only for k-mers shared with another haplotype or copy, so the change is on that other copy. `unresolved`: this sequence has too few k-mers of its own to tell |
| `confidence` | `high` or `low` (per class; see below) |
| `score` | class-specific strength: log-likelihood ratio in nats (S1, S3), reads (S2), unsupported k-mers (false joins), clipped reads (S4) |
| `channel` | S1 copy number, S2 missing k-mers, S3 parental hapmers, S4 long-read join audit |
| `summary` | one sentence a curator can read |
| `evidence` | `key=value` pairs with the numbers behind the call |

**`issues.bed`**: the same loci as BED9 coloured by class. Columns 1–4 are contig, start, end and
class, so the file can be given as it is to scorers that read a label in column 4.

| Class | Channel | What the evidence says | `low` confidence when |
|---|---|---|---|
| `COLLAPSE` | S1 | the reads hold more copies of this sequence than the release: a collapsed duplicate, or the other haplotype's copy is missing | a single depth factor explains the segment as well as a whole extra copy does (the issue then names what supports it, `support=`), or the segment has fewer than 5 assessed windows (`low_because=short`) |
| `EXPANSION` | S1 | the release holds one more copy than the reads support: a false duplication | as above (one copy fewer and half the depth cannot be told apart when only one k-mer class is present) |
| `UNSUPPORTED` | S1 | sequence the release holds once has no read support | as above |
| `BASE_ERROR` | S1 + S2, S1 + S4 | a run of up to k assembly k-mers without read support inside normally covered sequence, where the short reads carry another allele, or where long reads run through without carrying those k-mers: a wrong base or small error | fewer than 6 tracked k-mers, or fewer than 3 long reads running through |
| `MISSING_JUNCTION` | S2 | reads join this position to a distant one; no sequence in the release holds that junction | fewer than 6 reads or 5 tracked k-mers |
| `MISSING_HAPLOTYPE` | S2 | three or more alleles carried by the reads, less than 10 kb apart, are absent: one haplotype of the stretch is not represented | fewer than 5 alleles |
| `MISSING_ALLELE` | S2 | reads carry one allele here that no sequence in the release holds | fewer than 6 reads or 5 tracked k-mers |
| `MISSING_SEQUENCE` | S2 | reads run from here into sequence found nowhere in the release (reported with at least 6 reads and 5 tracked k-mers, or when its k-mers are in neither parent's reads; weaker ones stay in `missing_clusters.tsv`) | always |
| `PHASE_SWITCH` | S3 | a stretch carries the other parent's hapmers than the rest of its sequence, where both parents' reads cover it at their usual depth | fewer than 10 hapmer sites |
| `MISPLACED_PIECE` | S3 | a whole gap-delimited scaffold piece carries the other parent's hapmers | never |
| `UNSUPPORTED_JOIN` | S4 | no long read links the two sides of a run of N | unless another issue lies within 20 kb (a misplaced piece, a copy too many or too few, a missing junction): an unlinked scaffold gap alone is ordinary |
| `FALSE_JOIN` | S1 + S4 | a run of assembly k-mers that neither short nor long reads hold, and that no long read runs through (needs `--long-read-bam`) | fewer than 3 long reads clipped there and fewer than 10 over it |

Evidence keys of copy-number issues: `delta` (copies in the reads minus copies in the release),
`copies_by_multiplicity` (copies the reads give to k-mers the release holds 1, 2, ... times),
`own_specific` (the same for k-mers found only in this sequence's group), `fit_cn` and `fit_depth`
(the worst z-score of a k-mer class under the whole-copy and the depth-factor reading),
`other_copies` (where else the segment's k-mers occur: sequence:start-end:shared k-mers),
`unique_at_copies` (one-copy k-mers at those copies, and the copies the reads give them) with
`copies_llr` (their log-likelihood for full depth over the segment's depth ratio), and, for
low-confidence issues, `low_because` (`depth_factor_fits`, then with `support` = `own_kmers`,
`other_channel`, `long_and_deep` or `unlinked_join`; or `short`). One event
shows at every copy of its sequence: a false duplication gives an issue on the extra copy, on the
copy it repeats, and a `partner` on the other haplotype; `other_copies` ties them together.

## What could not be assessed

- **`cannot_assess.bed`**: runs of windows where the copy-number channel makes no statement: too
  few usable k-mers, more than half of the k-mers above the multiplicity cap (satellites and other
  high-copy repeats), extreme GC, or mostly N. Its size is reported in `summary.json`.
- **`depth_anomalies.tsv`**: segments whose depth departs from the release without fitting a whole
  copy: `DEPTH` (every k-mer class scaled by one factor: uneven sequencing depth, or sequence
  present in a fraction of cells) and `IRREGULAR` (no simple pattern). These are not issues. They
  are worth a look when the reads and the assembly come from different DNA.
- In `hapmer_sequences.tsv`, sequences with fewer than 20 hapmer sites are `not_assessable`, never
  "consistent".

## Per-channel tables

| File | Content |
|---|---|
| `cn_segments.tsv` | every copy-number segment with its kind (`CN`, `CN_OR_DEPTH`, `ABSENT`, `DEPTH`, `IRREGULAR`), best whole-copy change `delta`, best depth `factor`, copies implied by the reads per k-mer multiplicity (`classes`, as `multiplicity:copies`), the same for k-mers found only in this sequence's own group (`own_classes`), the other copies of the segment (`copies`), the one-copy k-mers at those copies (`copy_unique`, `copy_llr`) and whether they changed the kind (`resolved`) |
| `copies.hits.tsv`, `copies.regions.bed` | where the tracked k-mers of each candidate segment occur in the release (query = row of the segment table) |
| `unsupported_clusters.tsv` | every run of tracked assembly k-mers with no read support, with the support of its flanks and how many of its k-mers occur once in the release. A run is an issue (`BASE_ERROR`) only when the reads carry another allele there; a false join looks the same but so do about a thousand correct stretches of a truth genome that short reads do not cover |
| `missing_clusters.tsv` | clusters of solid read k-mers absent from the release: kind, locus, reads, tracked k-mers, median count, share found in each parent |
| `missing/` | the recruited reads (`missing.reads.fa`) and where their k-mers anchor (`missing.anchors.tsv`) |
| `hapmer_sequences.tsv`, `hapmer_pieces.tsv`, `hapmer_switches.tsv` | hapmer sites per sequence (with `windows_unreadable`: windows set aside because a parent's reads do not cover them at their usual depth), per gap-delimited piece, and the switch segments |
| `joins.tsv` | with `--long-read-bam`: per run of N, reads that link it, reads near it on each side, clips, verdict |
| `unsupported_points.tsv` | with `--long-read-bam`: per unsupported run with normally covered flanks (`run...`) and per copy-number loss on a haplotype's own k-mers (`cn...`): long reads that run through, long reads clipped there, the verdict by alignment (`spanned`, `broken`, `unclear`), the unread assembly k-mers, how many of them at least two long reads hold, and the verdict by content (`long_read_kmers`: `carried`, `not_carried`, `unclear`) |
| `cards/` | per issue, `ID.json` (the numbers) and `ID.png` (copies implied by the reads along the locus, hapmer sites, and where reads and release disagree in sequence) |

## Summary and provenance

- **`summary.json`**: depth per copy and the fitted noise model, assessable fraction, k-mer support
  of the assembly (with the consensus QV it implies), completeness of the read k-mers, counts of
  segments and issues by class, missing k-mers by copy class and by parent, timings.
- **`run.json`**: inputs with a checksum of their first megabyte, engine version, every parameter.
- **`model.json`**: depth per copy, the GC curve, extra-Poisson noise per k-mer multiplicity and the
  relative wander of depth.

## Engine files

`scan.track.bin` (one record per tracked k-mer position: gap to the previous position, multiplicity
per assembly group, count per read set), `scan.contigs.tsv`, `scan.gaps.bed`, `scan.gc100.u8`,
`scan.acgt100.u8`, `scan.novel.bin` (read k-mers absent from the release), `scan.spectrum.tsv`,
`scan.scan.json`, and `reads.*.cpkdb` (k-mer databases, reusable for any assembly of the same
sample). `copyproof.track.Scan` reads them.
