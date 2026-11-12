# HPRC release-2 assembly QC baselines (assembly review track)

Reproduces the measured baselines cited in the assembly review from public index tables in
github.com/human-pangenomics/hprc_intermediate_assembly (branch main, data_tables/).

Run: `scripts/hprc_r2_qc_baseline.sh <out_dir> [raw_base_url] [n_io_workers]`
(bash wrapper creates a self-cleaning `mktemp -d` scratch directory with `trap 'rm -rf' EXIT`;
python3 standard library only; ~70 MB of per-haplotype BED/TSV files are fetched over HTTPS from
s3-us-west-2.amazonaws.com/human-pangenomics and deleted on exit).

Inputs (index CSVs; md5 recorded in provenance.json at run time):
- data_tables/annotation/chrom_assignment/t2t_sequences_hprc_r2_v1.0.index.csv
- data_tables/annotation/chrom_assignment/gaps_hprc_r2_v1.0.index.csv
- data_tables/assembly_qc/flagger/flagger_hifi_v0.1.csv  (HMM-Flagger v1.2.0 conservative, HiFi)
- data_tables/assembly_qc/flagger/flagger_ont_v0.1.1.csv (HMM-Flagger v1.2.0 conservative, ONT)
- data_tables/assembly_qc/nucflag/nucflag_pre_release_v0.2.index.csv (NucFlag v0.3.3, HiFi)
- data_tables/sequencing_data/data_sample_level_summary_release2_v1.0.csv

Outputs: hprc_r2_qc_per_haplotype.tsv (462 rows), hprc_r2_t2t_by_chrom.tsv, provenance.json.
Derived in the review session: hprc_r2_qc_per_sample.tsv, hprc_r2_t2t_rate_by_chrom.tsv, hprc_r2_qc_summary.json.

Definitions: flagged bp = union of HMM-Flagger Err/Dup/Col intervals (NNN gap and Hap rows excluded);
HiFi-and-ONT bp = base-level intersection of the two unions on the same haplotype; NucFlag non-HET =
all NucFlag classes except HET (MISJOIN, COLLAPSE, COLLAPSE_VAR, COLLAPSE_OTHER). T2T levels as in the
HPRC t2t_sequences tables (contig = gapless, scaffold = T2T with >=1 gap). chrX/chrY denominators use
hprc_release2_sample_metadata.csv sex (115 male, 116 female among the 231 QC'd samples).
Caveats: 4 of 466 release-2 haplotypes (incl. HG002, Verkko) are absent from the QC indexes; base-level
overlaps are not window-expanded (HPRC2 reports 3.1% unreliable per haplotype using 500-kb windows and a
2-of-3 rule); NucFlag was designed for centromeres and its MISJOIN class includes coverage dips.
