# copyproof notebook

A running log of what was tried, what was seen and what changed, newest entry first. The guide proposes; this log records what the data decided.

## Entry template

```
### YYYY-MM-DD: short title
Question:  what this was meant to answer (link the GUIDE section or decision, e.g. D2)
Run:       commands or script; inputs (ids from resources.tsv); tool versions; where outputs live (outside git if large)
Seen:      the numbers that matter, with units and denominators
Decision:  what changes in the plan, or "no change"
Next:      the next step
```

## Open decisions

| # | Decision (GUIDE §9) | Status | Decided in |
|---|---|---|---|
| D1 | Read source for copy number | open | |
| D2 | k-mer size | open | |
| D3 | Coordinate frame | open | |
| D4 | Segmentation | open | |
| D5 | Negative-control denominator | open, settle before rung 3 | |
| D6 | Thresholds and score | open, freeze before rung 6 | |
| D7 | Non-trio assemblies | open | |
| D8 | Reuse or build | open | |
| D9 | Cost | open | |
| D10 | Change of course | open | |

## Baselines to reproduce

These values come from the pilot work in this repository. M0 and M1 should reproduce them before anything new is built on them.

| What | Value | Source |
|---|---|---|
| HG01943 read copy number in the duplicated segment (k-mers; CRAM index) | 3.13; 3.05 | report F4; `kmer_cn_track.py`, `crai_cn.py` |
| HG01943 assembly copy number in 10-kb bins | 87 bins at CN 3, 17 at CN 2, 3 at CN 4 | report F4 |
| HG01943 junction 31-mers in HPRC R2 | 0 of 28 | `regression_tests.md`, R1 |
| HG01943 piece B hapmers, maternal vs paternal | 5,638 vs 0 | R1 |
| HG01943 joins spanned by any read | 0 of 4 | R2 (T1); report F4 |
| HG01943 junction reads ending at the two ends of the maternal copy | 31 and 30 | report F4 |
| R2 targets hit: HMM-Flagger HiFi, conservative / full | 0/5 / 2/5 | R2 |
| R2 targets hit: HMM-Flagger ONT, conservative / full | 0/5 / 2/5 | R2 |
| R2 targets hit: NucFlag | 2/5 | R2 |
| NucFlag genome-wide MISJOIN calls, HG01943 paternal / maternal | 1,054 (about 361 per Gb) / 1,220 | R2 |
| HG01103 junctions in the trio assembly FASTA / seen via the published BAM | 16 of 16 / 5 of 16 | R3, R4 |
| HG01103 lost segments: paternal/maternal haplotagged depth | 0.42–0.46 (kept segments 0.95–0.97) | R5 |
| Case driver run time at 2 threads | HG01943 9 min 48 s; HG01103 8 min 2 s | `regression_tests.md` |
| Known low-confidence intervals in hg002v1.1 | 3,172 intervals, 38.8 Mb | resources.tsv: `hg002_known_issues_v1.1` |
| Samples with both parents in NYGC 30x | HPRC R2 121 of 223; HGSVC3 30 of 63 | `../data/sample_overlap.tsv` |

## Log

### 2026-11-13: guide written, resources checked

Question: what public data can drive each rung of the validation ladder (GUIDE §8)?

Run: listings of the HPRC data tables, the HGSVC3 release and working folders, and the HG002 Q100 folders. Every entry in resources.tsv was checked on this date.

Seen:

1. HPRC publishes, for every R2 haplotype, gap BEDs, segdup and censat annotation, chains to GRCh38 and CHM13, Flagger's read-to-assembly BAMs and high-clip bigWigs. S3–S5 can therefore start on HPRC R2 without realignment.
2. A known-issue BED exists for hg002v1.1 (3,172 intervals, 38.8 Mb). Rung 3 can separate known truth issues from new flags.
3. hg002v1.2 is public, with a chain to v1.1. Loci fixed between the versions are natural positives.
4. 24 HG002 assemblies have GQC outputs against v1.1, giving real, labelled errors for rung 4.
5. HGSVC3 ships unassigned and disconnected sequence in separate files. S1 needs to count it, or it will report false collapses.

Decision: none yet; these are inputs for M0–M2.

Next: M0. Re-run both case drivers and checksum the HG002 inputs.
