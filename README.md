# assemble_align

**Raising the state of the art in read alignment and genome assembly: a 2026 review, gap assessment and litmus pilots on public 1000 Genomes data.**

The main deliverable is the report: [`report/assemble_align_review.pdf`](report/assemble_align_review.pdf). Its markdown source is [`report/assemble_align_review.md`](report/assemble_align_review.md).

A one-page summary of the main findings is in [`docs/index.html`](docs/index.html), served at <https://jlanej.github.io/assemble_align/> once GitHub Pages is enabled for the `docs/` folder.

## Bottom line

Speed and contiguity are no longer the binding constraints. The remaining errors sit in segmental duplications, satellites and structurally variable loci, and they are increasingly silent: aligners give wrong placements their top mapping quality, and assemblies carry misplaced copies and cell-line rearrangements that pass coverage-based QC. The highest-impact opportunity is **calibrated placement** of reads and contigs, trained and tested against the HPRC release-2 and HGSVC3 haplotypes of 1kG samples.

The recommended order of work:

1. An open placement-truth benchmark on 1kG genomes (F1).
2. Calibrated MAPQ with paralog posteriors (F2).
3. Placement-aware assembly QC (F4).
4. A minimap2 sketching patch (F6).
5. Then structural population-aware references (F3), assembler repairs (F4, F5) and population-scale ONT assembly (F7).

## Layout

| Path | Contents |
|---|---|
| `docs/` | Landing page for copyproof (`index.html`), with six figures drawn from the case tables (`assets/fig_cp_*.svg`); the one-page summary of the review (`review.html`, five figures in `assets/`) |
| `report/` | Report PDF and markdown, figures, `ranked_proposals.tsv`, and the build script (`build_pdf.sh`: pandoc → HTML → weasyprint) |
| `references/` | `references.bib` (369 cited works) and `references.json` (452 verified records, with a `cited_in_report` flag); `tracks/` holds each review track's own bibliography |
| `reviews/` | Full literature reviews: alignment, assembly, transferable ideas (`ideas_matrix.tsv`); verified tool versions as of 2026-11-11; the HPRC release-2 QC baseline script and tables |
| `data/` | Catalog of 61 public 1kG-linked datasets, a 2,122-sample overlap matrix and its builder, technology error profiles, the open benchmark protocol, and remote-access recipes |
| `litmus/L1_seed_sampling/` | Low-density sketching in minimap2: numba port of `mm_sketch`, ten sampling schemes, and the `--mod-oc` patch (`minimap2_sampling.patch`) |
| `litmus/L2L3_reference_headroom/` | NA19240 chr20 real reads (NYGC Illumina, HPRC HiFi/ONT) against five references, including the oracle and personalized references |
| `litmus/L5_mapq_calibration/` | MAPQ calibration with truth lifted through the HPRC impg alignments; a learned recalibrator with an apply-to-BAM script |
| `litmus/L4_assembly/` | Assembly pilot on real reads: chr20 coverage and read-type titration with hifiasm 0.25.0 (HG01258; HiFi, R10 ONT, HiFi + R10 as `--ul`) and four read-recruitment rules at 13 hard loci (NA19240); scripts 01-14, tables, figures |
| `factcheck/` | Independent fact-check of the draft: report, 52 findings, how each was resolved, and spot-check scripts |
| `litmus/cases_HG01943_HG01103/` | Re-test of the two flagged 1kG cases, five scored regression tests, drivers and scorer |
| `copyproof/` | copyproof, a placement-aware assembly auditor (the detection half of F4): a working prototype (k-mer engine in Rust, Python package, tests), the implementation guide it started from, the method as built, a running notebook, and one record per experiment on public data |

## Headline litmus results

- **L5, MAPQ calibration (simulated NA19240 reads).** Inside segmental duplications, MAPQ 60 is wrong for 1.23-2.77% of reads on GRCh38 and 0.49-1.42% on CHM13. (Long reads are scored against ≥ 50%-lifted truth; under a strict ≥ 90%-colinear truth their in-SD error is 0.12-0.33%.) A chromosome-held-out recalibrator cuts in-SD calibration error by 59-93%. Ranked by its score, with the threshold set on held-out labels, it keeps 28-58% of correct in-SD reads at ≤ 0.1% error, a level no reported-MAPQ threshold reaches; label-free thresholds at a 0.1% target realise 0.06-0.40% in-SD error, against 2.8-10.5% for reported MAPQ.
- **L1, sketching (simulated reads).** A strand-symmetric open-closed mod-minimizer cuts minimap2 index entries by 13.1-26.1% without reducing the share of unique-sequence long reads with ≥ 3 true-origin seeds (−0.11 points for Illumina). The 241-line patch gives 14-22% smaller indexes and 12-30% less mapping CPU on chr20 at equal or lower wrong-mapping rates at MAPQ ≥ 1, at 1.3-1.9x slower index builds.
- **L2/L3, reference headroom (real reads, chr20).** Small-variant personalization captures 88-93% of the Illumina unique-sequence gap (own calls 88%, 1kGP phased 93%) but only 9-11% of the SD reference-bias excess. The remaining headroom is structural.
- **Cases.** HG01943's de novo 1.09-Mb tandem duplication is missing from the HPRC release-2 assembly (0/28 junction 31-mers), and the conservative published QC tracks flag 0/5 targets. HG01103's sub-clonal chromothripsis is present in the trio assembly (16/16 junctions), but long contigs follow the minority clone at 11/16 forks.
- **L4, assembly (real reads; provisional, 6 of 15 planned runs).** On HG01258 chr20, hifiasm reproduces 99.6% of the HPRC truth at ≥ 99.9% identity from 36.7x HiFi and 98.3% from 35.8x R10 nanopore reads alone, but ONT-only falls to 73.5% at 16x (HiFi 97.6%), below what the 16.9x Vienna cohort could support. Standard-length R10 reads as `--ul` cut HiFi contigs from 16 to 5. Population k-mers restore oracle-level recruitment at HLA (21 → 24/26 resolved haplotype-loci) but do not beat alignment-window recruitment; no rule resolves SMN.

## Reproducing

Each `litmus/*/` folder has a `README.md` (the full write-up), scripts and small result tables. Raw data are never stored in the repository. Scripts stream region slices over HTTPS from public buckets: the 1000 Genomes EBI FTP and S3, `1000g-ont`, `human-pangenomics` (HPRC), and HGSVC3. They write intermediates into `mktemp -d` scratch directories that are removed on exit (`trap ... EXIT`). Tools used: minimap2 2.31, samtools/bcftools 1.2x, bwa 0.7.19, hifiasm 0.25.0, yak, meryl, numba, scikit-learn; versions are recorded in each pilot's provenance tables. The pilots ran on a 10-core, 64 GB laptop; single-part whole-genome minimap2 indexes need about 13-16 GB of RAM.

To rebuild the PDF: `cd report && ./build_pdf.sh assemble_align_review.md assemble_align_review.pdf` (needs pandoc ≥ 3 and weasyprint).

## Caveats

- L1 and L5 use simulated reads; the real-read acceptance tests are specified in the report (F1, F2, F6).
- L2/L3 covers one sample and one chromosome.
- L4 is provisional: a partial titration grid on one sample (HG01258, chr20), hard loci on one sample (NA19240), hifiasm `-f0`, and chromosome-sliced Illumina k-mers for QV and trio phasing.
- The case studies use region slices.
- Several 2025-2026 sources are preprints.
- The proposal scores are structured judgments.
- The re-test corrects two details of the original case repositories (HG01943: all 31 junction reads carry a 2-bp `CT` insertion; HG01103: the trio paternal assembly holds all 16 junctions, and the published alignment BAM omits contigs and contig ends); see `litmus/cases_HG01943_HG01103/`.
- The independent fact-check's 52 findings are applied; `factcheck/factcheck_resolution.tsv` records how each was resolved.

## License

- **Code** (all `*.py`, `*.sh`, `*.c` and `*.lua` files, `report/style.css` and `report/build_pdf.sh`, and `litmus/L1_seed_sampling/minimap2_sampling.patch`): MIT License, see [`LICENSE`](LICENSE). The patch modifies minimap2, which is MIT-licensed by its authors; the minimap2 portions keep their original copyright notice.
- **Report, reviews, figures, write-ups and result tables**: Creative Commons Attribution 4.0 International, see [`LICENSE-CC-BY-4.0.md`](LICENSE-CC-BY-4.0.md).
- **Third-party data.** Tables derived from public 1000 Genomes, HPRC, HGSVC, GIAB and 1kGP-LRSC data remain subject to those projects' data-use and citation terms. Cite the original data releases (listed in `references/`) when reusing them.

The report was prepared on 2026-11-11/12 from public data only.
