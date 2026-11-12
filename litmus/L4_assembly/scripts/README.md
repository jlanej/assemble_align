# L4 assembly litmus pilot: scripts

The pipeline runs on public 1000 Genomes / HPRC data. Every script that writes intermediates uses `mktemp -d "${TMPDIR:-$PWD}/scratch.XXXXXX"` with a `trap ... EXIT` cleanup and takes its inputs (URLs or paths) and outputs as arguments. Tested on macOS arm64 with hifiasm 0.25.0-r726, minimap2 2.31-r1302, samtools/htslib 1.25, yak 0.1-r69, KMC 3.2.4, seqkit 2.14.0 and Python 3.11 (pandas, numpy, matplotlib). Thread budget: 3.

| Step | Script | What it does |
|---|---|---|
| 1 | `01_fetch_longread_slice.sh BAM_URL REGION OUT` | primary alignments in a region of a remote indexed BAM -> FASTQ (chr20 HiFi/ONT recruitment) |
| 2 | `02_fetch_truth_seq.sh FA_GZ_URL OUT.fa SEQ[:s-e]...` | remote `samtools faidx` of HPRC release-2 haplotypes (truth chr20) |
| 3 | `03_fetch_illumina_slice.sh CRAM CRAI REF REGION OUT` | NYGC 1kGP CRAM slice -> FASTQ (child and parents) |
| 4 | `04_locus_haplotypes.py loci.tsv assemblies.tsv OUTDIR [--truth S1,S2] [--only LOCUS]` | project GRCh38 locus windows onto HPRC r2 haplotypes via minigraph-cactus chains (overlap or unique-anchor mode) and fetch sequences |
| 5 | `05_nested_subsample.py reads.fq.gz G OUT 8,16,full` | nested coverage subsamples (hash-ordered) |
| 6 | `06_run_hifiasm.sh hifi/hifi_ul/ont OUTDIR THREADS READS [UL]` | hifiasm 0.25.0 with `-f0`; exports hap1/hap2/primary FASTA and the run's time/RSS line |
| 7 | `07_eval_asm.py --asm DIR --pat --mat --child --patyak --matyak --label --out` | cov999, NG50, NGA50, yak QV, yak trioeval, truth-window phasing |
| 8 | `08_fetch_pool_bam.sh BAM_URL pool.bed OUT` | candidate pool for locus recruitment (windows +/- 1-2 Mb plus unmapped reads) |
| 9 | `09_loci_pipeline.py --sample ...` | REF_ALN / REF_KMER / POP_KMER / ORACLE recruitment -> hifiasm -> spanned/resolved vs truth |
| 10 | `10_yak_dbs.sh OUT THREADS child.fq father.fq mother.fq` | yak k=31 databases |
| 11 | `11_titration_driver.sh S MODE READS G COVS WORKDIR THREADS [UL]` | subsample -> assemble -> evaluate loop |
| 12 | `12_tables_figures.py WORKDIR OUTDIR` | L4_titration.tsv, L4_loci.tsv, two PNG figures |
| 13 | `13_master_queue.sh WORKDIR SCRIPTS_DIR [THREADS]` | sequential, re-entrant end-to-end queue; skip markers (`hap1.fa` + `eval.jsonl` + `SKIPPED`) drop titration points |
| 14 | `14_build_report.py OUTDIR methods.md fig1 fig2` | composes the markdown results with computed numbers |

Inputs used here: see `L4_provenance.tsv` (URLs and md5) and `L4_loci_definitions.tsv`. The truth is HPRC release-2 trio-phased haplotypes (`*_hprc_r2_v1.0.1.fa.gz`). chr20 contigs come from `annotation/chrom_assignment/*chromAlias.txt`, and locus haplotypes from `annotation/chains/minigraph-cactus/*_vs_GRCh38.chain.gz`.
