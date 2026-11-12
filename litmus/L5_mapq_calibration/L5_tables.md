### Table 1. Reported MAPQ versus empirical error

| Run | Evaluable reads | MAPQ 0: n (error) | MAPQ 60: n, wrong | MAPQ 60 error, outside SD | MAPQ 60 in SD: wrong / n | MAPQ 60 error, in SD (Phred) | MAPQ≥30 error, in SD |
|---|---|---|---|---|---|---|---|
| GRCh38 · Illumina · bwa mem | 1,163,145 | 33,700 (75.6%) | 1,084,363, 605 | 0.00019 | 403 / 32,777 | 1.23% (19) | 5.44% |
| GRCh38 · Illumina · minimap2 sr | 1,162,960 | 30,460 (73.8%) | 1,083,298, 860 | 0.00038 | 458 / 34,160 | 1.34% (19) | 3.65% |
| GRCh38 · HiFi · minimap2 map-hifi | 37,114 | 100 (57.4%) | 33,877, 140 | 9.5e-05 | 130 / 8,261 | 1.57% (18) | 2.58% |
| GRCh38 · ONT · minimap2 lr:hq | 21,651 | 25 (48.0%) | 20,357, 141 | 4.9e-05 | 136 / 4,906 | 2.77% (16) | 3.55% |
| GRCh38 · ONT · minimap2 map-ont | 21,651 | 26 (53.8%) | 20,293, 136 | 4.9e-05 | 131 / 4,842 | 2.71% (16) | 3.44% |
| CHM13 · Illumina · bwa mem | 1,181,658 | 53,200 (80.9%) | 1,085,570, 276 | 0.00011 | 158 / 32,277 | 0.49% (23) | 2.48% |
| CHM13 · Illumina · minimap2 sr | 1,180,566 | 49,322 (80.2%) | 1,083,671, 284 | 0.00011 | 166 / 33,297 | 0.50% (23) | 1.79% |
| CHM13 · HiFi · minimap2 map-hifi | 38,813 | 166 (73.2%) | 33,962, 92 | 7.2e-05 | 86 / 8,326 | 1.03% (20) | 1.73% |
| CHM13 · ONT · minimap2 lr:hq | 22,741 | 67 (91.2%) | 20,458, 75 | 4.8e-05 | 70 / 4,933 | 1.42% (18) | 2.37% |
| CHM13 · ONT · minimap2 map-ont | 22,741 | 76 (72.2%) | 20,373, 69 | 4.8e-05 | 64 / 4,843 | 1.32% (19) | 2.23% |

### Table 2. Held-out recalibration inside segmental duplications

| Run | ECE reported | ECE MAPQ-only | ECE GBM | ECE GBM+iso | ECE change (GBM+iso) | Brier reported → GBM+iso | AUPRC(wrong) reported → GBM | Correct kept at ≤0.1% error: reported / GBM / GBM+iso |
|---|---|---|---|---|---|---|---|---|
| GRCh38 · Illumina · bwa mem | 0.183 | 0.055 | 0.015 | 0.013 | -93% | 0.184 → 0.104 | 0.57 → 0.72 | 0% / 45% / 43% |
| GRCh38 · Illumina · minimap2 sr | 0.150 | 0.058 | 0.019 | 0.019 | -87% | 0.164 → 0.105 | 0.57 → 0.73 | 0% / 40% / 36% |
| GRCh38 · HiFi · minimap2 map-hifi | 0.065 | 0.035 | 0.041 | 0.011 | -83% | 0.076 → 0.064 | 0.41 → 0.42 | 0% / 56% / 0% |
| GRCh38 · ONT · minimap2 lr:hq | 0.065 | 0.041 | 0.039 | 0.015 | -77% | 0.073 → 0.063 | 0.37 → 0.41 | 0% / 28% / 0% |
| GRCh38 · ONT · minimap2 map-ont | 0.064 | 0.040 | 0.042 | 0.016 | -76% | 0.073 → 0.065 | 0.36 → 0.39 | 0% / 42% / 0% |
| CHM13 · Illumina · bwa mem | 0.154 | 0.067 | 0.019 | 0.014 | -91% | 0.157 → 0.098 | 0.59 → 0.71 | 0% / 46% / 46% |
| CHM13 · Illumina · minimap2 sr | 0.145 | 0.070 | 0.025 | 0.022 | -85% | 0.152 → 0.098 | 0.59 → 0.73 | 0% / 45% / 35% |
| CHM13 · HiFi · minimap2 map-hifi | 0.048 | 0.035 | 0.018 | 0.013 | -73% | 0.064 → 0.057 | 0.55 → 0.63 | 0% / 58% / 0% |
| CHM13 · ONT · minimap2 lr:hq | 0.049 | 0.027 | 0.024 | 0.020 | -60% | 0.061 → 0.057 | 0.52 → 0.58 | 0% / 38% / 0% |
| CHM13 · ONT · minimap2 map-ont | 0.050 | 0.030 | 0.021 | 0.021 | -59% | 0.062 → 0.057 | 0.52 → 0.58 | 0% / 40% / 0% |

### Table 3. Selecting reads with the score alone (thresholds never see truth labels)

| Run | in SD: MAPQ≥30 error / kept, reported | in SD: MAPQ≥30 error / kept, GBM+iso | in SD: realised error at 0.1% target, reported / GBM / GBM+iso | genome-wide: MAPQ≥30 error / kept, reported | genome-wide: MAPQ≥30 error / kept, GBM+iso | genome-wide: realised error at 0.1% target, reported / GBM / GBM+iso |
|---|---|---|---|---|---|---|
| GRCh38 · Illumina · bwa mem | 5.44% / 76% | 0.05% / 35% | 10.46% / 0.20% / 0.06% | 0.37% / 98% | 0.00% / 95% | 1.16% / 0.15% / 0.13% |
| GRCh38 · Illumina · minimap2 sr | 3.65% / 72% | 0.10% / 36% | 5.71% / 0.17% / 0.20% | 0.23% / 98% | 0.01% / 95% | 0.62% / 0.16% / 0.16% |
| GRCh38 · HiFi · minimap2 map-hifi | 2.58% / 88% | 0.11% / 53% | 4.13% / 0.51% / 0.22% | 0.14% / 99% | 0.00% / 97% | 0.39% / 0.25% / 0.10% |
| GRCh38 · ONT · minimap2 lr:hq | 3.55% / 91% | 0.24% / 44% | 5.00% / 0.80% / 0.40% | 0.18% / 99% | 0.01% / 95% | 0.38% / 0.25% / 0.13% |
| GRCh38 · ONT · minimap2 map-ont | 3.44% / 90% | 0.28% / 39% | 4.67% / 0.75% / 0.33% | 0.18% / 99% | 0.01% / 97% | 0.39% / 0.26% / 0.11% |
| CHM13 · Illumina · bwa mem | 2.48% / 74% | 0.04% / 30% | 6.15% / 0.19% / 0.09% | 0.18% / 98% | 0.00% / 94% | 0.76% / 0.15% / 0.10% |
| CHM13 · Illumina · minimap2 sr | 1.79% / 70% | 0.08% / 35% | 3.35% / 0.19% / 0.13% | 0.11% / 98% | 0.00% / 95% | 0.44% / 0.14% / 0.13% |
| CHM13 · HiFi · minimap2 map-hifi | 1.73% / 87% | 0.17% / 40% | 2.79% / 0.10% / 0.14% | 0.12% / 99% | 0.00% / 96% | 0.46% / 0.17% / 0.07% |
| CHM13 · ONT · minimap2 lr:hq | 2.37% / 89% | 0.25% / 41% | 3.47% / 0.46% / 0.23% | 0.15% / 99% | 0.01% / 95% | 0.39% / 0.21% / 0.14% |
| CHM13 · ONT · minimap2 map-ont | 2.23% / 88% | 0.34% / 50% | 3.21% / 0.30% / 0.34% | 0.14% / 99% | 0.02% / 97% | 0.44% / 0.23% / 0.13% |

### Table 4. Why MAPQ-60 placements are wrong (each read realigned to its truth and reported loci ±2 kb)

| Run | Wrong at MAPQ 60 | In SD | Read scores better at reported locus | Tie | Truth locus better (search error) | Read does not align at truth locus | Median truth-SD identity |
|---|---|---|---|---|---|---|---|
| GRCh38 · Illumina · bwa mem | 605 | 67% | 501 (83%) | 41 | 3 | 59 | 0.984 |
| GRCh38 · Illumina · minimap2 sr | 860 | 53% | 659 (77%) | 93 | 19 | 89 | 0.979 |
| GRCh38 · HiFi · minimap2 map-hifi | 140 | 93% | 120 (86%) | 1 | 19 | 0 | 0.995 |
| GRCh38 · ONT · minimap2 lr:hq | 141 | 96% | 134 (95%) | 0 | 7 | 0 | 0.997 |
| GRCh38 · ONT · minimap2 map-ont | 136 | 96% | 128 (94%) | 0 | 8 | 0 | 0.997 |
| CHM13 · Illumina · bwa mem | 276 | 57% | 232 (84%) | 34 | 0 | 10 | 0.993 |
| CHM13 · Illumina · minimap2 sr | 284 | 58% | 201 (71%) | 70 | 4 | 9 | 0.987 |
| CHM13 · HiFi · minimap2 map-hifi | 92 | 94% | 85 (92%) | 0 | 7 | 0 | 0.995 |
| CHM13 · ONT · minimap2 lr:hq | 75 | 93% | 71 (95%) | 0 | 4 | 0 | 0.995 |
| CHM13 · ONT · minimap2 map-ont | 69 | 93% | 66 (96%) | 0 | 3 | 0 | 0.995 |

### Table 5. Long reads under a stricter truth (≥90% of the read colinear with its orthologous locus)

| Run | MAPQ 60 in-SD error: ≥50% / ≥90% colinear truth | ECE reported / MAPQ-only / GBM / GBM+iso (strict) | Correct kept at ≤0.1%: reported / GBM / GBM+iso (strict) | Realised error at 0.1% target: MAPQ-only / GBM+iso (strict) |
|---|---|---|---|---|
| GRCh38 · HiFi · minimap2 map-hifi | 1.57% / 0.33% | 0.035 / 0.011 / 0.024 / 0.005 | 0% / 71% / 0% | 0.41% / 0.15% |
| GRCh38 · ONT · minimap2 lr:hq | 2.77% / 0.12% | 0.024 / 0.009 / 0.019 / 0.009 | 0% / 71% / 0% | 0.19% / 0.29% |
| GRCh38 · ONT · minimap2 map-ont | 2.71% / 0.15% | 0.023 / 0.007 / 0.014 / 0.009 | 0% / 74% / 0% | 0.19% / 0.28% |
| CHM13 · HiFi · minimap2 map-hifi | 1.03% / 0.22% | 0.031 / 0.024 / 0.023 / 0.013 | 0% / 67% / 0% | 0.25% / 0.53% |
| CHM13 · ONT · minimap2 lr:hq | 1.42% / 0.23% | 0.022 / 0.016 / 0.020 / 0.011 | 0% / 70% / 0% | 0.27% / 0.25% |
| CHM13 · ONT · minimap2 map-ont | 1.32% / 0.14% | 0.022 / 0.019 / 0.023 / 0.011 | 0% / 74% / 0% | 0.16% / 0.45% |
