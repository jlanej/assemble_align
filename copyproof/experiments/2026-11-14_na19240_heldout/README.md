# Rung 6: NA19240, the held-out test

**State of this folder: the pre-registration.** This text was written at 10:22 on 2026-11-14
(CST). At that time the reads of the child and of the mother had been counted, the father's were
being counted, and no run of copyproof on NA19240 had started. The record of the test replaces
this text when the runs are done.

**Question.** With every threshold fixed beforehand, do the results hold on a second person
(GUIDE §8 rung 6)?

**What is fixed, and where.**

- `NOTEBOOK.md`, entry p: the sample and the release under test, the evidence (two other
  assemblies of the same person as a proxy for a truth genome), the three scores with their
  denominators, the rule for the verdict, what is reported without a prediction, and one run that
  was added before any read was counted. Each addition to the entry carries its time.
- `FROZEN_COMMIT`: the commit whose code and defaults are under test. The drivers refuse to run
  when `src/`, `engine/src` or `scripts/run_sample.py` differ from it.
- `predictions.tsv`: the seven counts that NA19240 is held against, from HG002 (written by
  `tables.py` after `calibrate_hg002.sh`).
- `prereg.sha256`: checksums of the scoring script, the drivers, the manifest, the predictions and
  each wording of entry p, with the time each line was added.

**What has been run.** The calibration on HG002, where the truth is known (`calibrate_hg002.sh`;
tables `hg002_*`, `score1_counts.tsv`, `alt_settings_hg002.tsv`). The downloads of
`scripts/manifests/05_na19240_heldout.tsv`, each file checked against its published MD5
(`inputs.md5`, `check_inputs.sh`). The k-mer counts of the reads and of the assemblies that serve
as label sets (`kmer_databases.tsv`). The alignment of the two haplotypes that the planted errors
need (`plant.sh`, stages `align` and `chain`).

**What has not been run.** Everything that produces a result on NA19240: the scored run and the
run with long reads (`scripts/run_na19240_heldout.sh`), the planted errors (`plant.sh`), the
hapmer audit on the HGSVC3 Verkko assembly (`hapmers_against_hprc.py`) and the alternative
settings (`alt_settings.py`).

**Checked on a stand-in.** `tables.py` writes the verdicts. Its half that reads the runs of
NA19240 was run once on the HG002 calibration run under NA19240's file names, with made-up QC
files, before any run of NA19240 existed: it returned the seven predicted counts, and "holds" for
each. Two faults in reading files were corrected on that occasion (the BED reader cut a line at
the first letter t, which the names in NA19240's files happen not to hold; an empty table of the
hapmer audit stopped the script), and the names of the label sets are now read from the score
file. The checksum of the corrected script is in `prereg.sha256` with the time.

**Not opened.** HPRC's Flagger and NucFlag calls for NA19240 are used to place the planted errors
and are compared with the issues afterwards. Their first three lines were read to learn the
format; beyond that they had been neither looked at nor used when this was written.
