# L1: near-optimal k-mer sampling for minimap2 — index size, anchor load and seed sensitivity

*Litmus pilot L1, 2026-11-11. Data: CHM13v2.0 and GRCh38 chr20, NA19240 HPRC release-2 haplotypes. All reads simulated with per-base truth. minimap2 v2.31 (r1302).*

**Bottom line.** A strand-symmetric open-closed mod-minimizer ([Groot Koerkamp et al. 2025](https://doi.org/10.1186/s13015-025-00270-0)), run at minimap2 2.31's unchanged (k, w), indexes 13–26 % fewer k-mers than minimap2's random minimizer at every preset tested:

| Preset | Index entries vs minimap2 |
|---|---|
| sr | −26.1 % |
| map-ont | −24.6 % |
| map-hifi / asm5 | −18.8 % |
| asm20 | −24.6 % |
| lr:hqae | −13.1 % |

GRCh38 chr20 reproduces each reduction to within 0.2 pp. Anchors per read fall by 14–23 % (8.5 % at lr:hqae). Every unique-sequence long read or haplotype tile keeps ≥ 3 true-origin seeds (Δ 0.00 pp); for Illumina reads the change is −0.11 pp. Unique-sequence reads are ≥ 94 % of NA19240-derived reads.

The gain is efficiency, not sensitivity: at equal density, long-read true-seed yield changes by ≤ 1.5 %.

A 241-line patch adds `--mod-oc` and persists it with a new index magic. It builds on macOS arm64. Against stock 2.31 on chr20 it gives 14–22 % smaller indexes and 12–30 % less mapping CPU, at equal or lower wrong-mapping rates (paftools mapeval). On the whole CHM13v2.0 genome the .mmi is 14.2 % (map-hifi) and 22.1 % (map-ont) smaller.

Two side findings:
1. minimap2's `hash64` is a weak random order for k ≤ 15. At the default map-ont preset it adds 2.8 % index entries, fixable with an invertible width-adapted mixer.
2. Naive canonicalization erases most of the gain of orders optimised on one strand (GreedyMini, double decycling). Open-closed schemes stay lossless when k − s is even.

The residual risk sits in repeat-dominated reads (0.6–5.8 %), where low-density schemes keep fewer sub-cutoff seeds. The seed metrics do not model minimap2's occ_dist rescue, but the aligner-level tests showed no accuracy loss.

![Density factor and index reduction per preset (a), true seeds (b) and longest seedless gap (c) vs index-entry change, and equal-density curves for NA19240-derived Illumina, ONT and HiFi reads (d–f). Panels b–c: NA19240 reads/tiles, one point per scheme × preset; c uses the unique-sequence stratum. Dashed line in d–f: minimap2's preset density.](L1_density_sensitivity.png)

## Litmus statement

**L1-A, low-density sampling in minimap2 (headline).** Replace minimap2 2.31's random minimizer with the strand-symmetric open-closed mod-minimizer (`--mod-oc`), with every preset's (k, w) unchanged.

*Measured baseline (CHM13 chr20, stock 2.31).*
* Index entries: sr 11,069,512; map-ont 12,452,943; map-hifi/asm5 6,652,482; asm20 12,050,520; lr:hqae 2,569,321. Density factors 2.006–2.069.
* .mmi size: sr 197.1 MB; map-ont 210.1 MB; map-hifi 131.0 MB.
* Mapping CPU, 1 thread: 5.23 s (79,281 Illumina reads), 37.7 s (1,574 ONT reads), 14.9 s (2,942 HiFi reads).
* Wrong mappings at MAPQ ≥ 1: 0.319 % / 0.259 % / 0.345 %.
* Whole-genome CHM13v2.0 map-hifi: .mmi 5.99 GB, build peak RSS 7.70 GB, wrong mappings at MAPQ ≥ 1 0.482 %.

*Measured effect of the prototype.*
* Index entries −26.1 / −24.6 / −18.8 / −24.6 / −13.1 % (GRCh38 chr20 within 0.2 pp).
* Anchors per read −8.5 % (lr:hqae tiles) and −13.8 to −22.6 % at the other presets.
* No loss in the fraction of unique-sequence reads with ≥ 3 true seeds for long reads (0.00 pp); −0.11 pp for Illumina.
* .mmi −14 to −22 %; mapping CPU −12 to −30 %; wrong mappings at MAPQ ≥ 1 equal or lower.
* Whole-genome map-hifi: .mmi −14.2 % and build RSS −33 %, but build CPU 1.9×.

*Target and pass/fail for the next, whole-genome test.* Use real HPRC/1kGP reads plus simulated reads with liftover truth, ≥ 10× coverage per technology. The prototype PASSES only if all of the following hold:
1. Index entries are reduced by ≥ 15 % at every preset (≥ 20 % where k > w).
2. Mapping CPU per Gb of reads is reduced by ≥ 10 %.
3. The wrong-mapping rate at MAPQ ≥ 1 and at MAPQ ≥ 10 rises by no more than +0.05 pp absolute.
4. Unmapped reads rise by no more than +0.1 pp overall and by no more than +1 pp in CenSat/SegDup strata.

Any violation is a FAIL.

**L1-B, hash quality (cheap fix).** Baseline: minimap2-exact density at w10/k15 is 0.18655 on i.i.d. sequence, 2.6 % above 2/(w+1). On CHM13 chr20 it is 0.18808, against 0.18277 for a strong random order.

Target: a width-adapted invertible 2k-bit mixer (`hash_w`) reaches ≤ 1.003 × 2/(w+1) on random sequence. It must cut k = 15 preset entries by ≥ 2 % on chr20, leave k ≥ 16 indexes unchanged (keep hash64 for 2k ≥ 32), and leave true seeds per read within ±0.5 % of the density-proportional expectation. Measured on random sequence: 0.18180 (PASS on the density criterion); the chr20 check is in L1_results.tsv (row hashw).

**L1-C, equal-density sensitivity for short reads.** Baseline: Illumina NA19240 reads, minimap2 sr, 17.7 true seeds per read and a 35-bp p95 longest gap.

Closed syncmers at matched density (+0.8 % entries) gave +9.7 % [9.5, 9.9] true seeds and a 31-bp p95 gap. Target for a mapping test: ≥ +0.1 pp correctly mapped reads in difficult regions (GIAB stratifications, SegDups) at ≤ +2 % index entries.


## Background

Minimap2 indexes a reference by its (w,k)-minimizers. Each window of w consecutive k-mers contributes the k-mer with the smallest hash, giving a density near 2/(w+1) ([Roberts et al. 2004](https://doi.org/10.1093/bioinformatics/bth408); [Schleimer et al. 2003](https://doi.org/10.1145/872757.872770); [Li 2018](https://doi.org/10.1093/bioinformatics/bty191)). Density sets index size and the number of seed lookups and anchors per read. The question is how much of it is unnecessary. Theory says a lot is. Forward schemes have a near-tight lower bound ([Kille et al. 2024, preprint](https://doi.org/10.1101/2024.09.06.611668), Theorem 2); the simple expression ⌈(w+k)/w⌉/(w+k) equals it only when k ≡ 1 (mod w) and otherwise understates it (at w = k = 19 the bound is about 3/39 ≈ 0.077), and density factors well below the random minimizer's 2.0 have long been known to be achievable ([Marçais et al. 2017](https://doi.org/10.1093/bioinformatics/btx235); [Marçais et al. 2018](https://doi.org/10.1093/bioinformatics/bty258)).

Practical schemes fall into two regimes.

- **Small k (k ≤ w).** These schemes first sample a context-free k-mer subset, then break ties by a random order. Miniception prefers closed syncmers ([Zheng et al. 2020](https://doi.org/10.1093/bioinformatics/btaa472)). Decycling-set orders use the Mykkeltveit embedding ([Pellow et al. 2023](https://doi.org/10.1101/gr.277644.123)). The open-closed minimizer prefers open syncmers and then closed syncmers ([Groot Koerkamp et al. 2025](https://doi.org/10.1186/s13015-025-00270-0)). Greedy and exact searches for low-density orders also belong here ([Golan et al. 2025](https://doi.org/10.1093/bioinformatics/btaf251); [Shur et al. 2026](https://doi.org/10.1007/978-3-032-42866-0_20)).
- **Large k (k > w).** The mod-minimizer selects the smallest t-mer, t = r + ((k−r) mod w), and samples the k-mer at its offset mod w. Its density tends to 1/w ([Groot Koerkamp & Pibiri 2024, preprint of WABI 2024](https://doi.org/10.1101/2024.05.25.595898)). Mod-sampling over an open-closed t-mer order, the open-closed mod-minimizer, improves on it at every k ([Groot Koerkamp et al. 2025](https://doi.org/10.1186/s13015-025-00270-0)).

Most of these papers measure density on random strings with forward-only schemes.

Sensitivity is a separate axis from density. Syncmers are context-free, so an unchanged k-mer keeps its sampling decision regardless of nearby mutations. That gives them higher conservation than minimizers at equal density ([Edgar 2021](https://doi.org/10.7717/peerj.10805); [Shaw & Yu 2021](https://doi.org/10.1093/bioinformatics/btab790); [Dutta et al. 2022](https://doi.org/10.1371/journal.pcbi.1010638)). Strobealign builds its seeds from syncmers ([Sahlin 2022](https://doi.org/10.1186/s13059-022-02831-7)). Repeat-aware down-weighting of k-mers is a further, orthogonal option ([Jain et al. 2020](https://doi.org/10.1093/bioinformatics/btaa435); [Jain et al. 2022](https://doi.org/10.1038/s41592-022-01457-8)).

Minimap2 2.31 still uses random minimizers. An open draft pull request (lh3/minimap2 #1361, 2026-08-25) adds a sketch-time mod-minimizer. It reports −19.5 % sampled positions on CHM13 chr21 at map-ont and non-inferior mapping accuracy. It is sketch-time only and refuses `-d` and prebuilt `.mmi` indexes, because its author found that a persisted mod-minimizer index without a new magic number is mapped silently wrong by old binaries; it gives no gain when t = k, i.e. at map-hifi, lr:hq, asm5 and asm10.

The gap this pilot measures has four parts:

1. Which scheme minimap2 should use at each preset's (k, w).
2. How much density survives the strand-symmetric (canonical) construction minimap2 requires.
3. Whether density savings cost true-origin seed matches under realistic errors and real haplotype divergence.
4. Whether a persisted index can be made safe.

## Methods

### minimap2 v2.31 sketching (read from source, tag v2.31 = r1302, tarball md5 ca5ddc443d95f040637fb7bc3ff17bb1)
* `mm_sketch()` (sketch.c): 2-bit codes A0 C1 G2 T3; forward and reverse-complement k-mer codes are maintained; the canonical strand is z = [fwd > rc]; palindromic k-mers (fwd == rc, only possible for even k) are skipped without taking a window slot; an ambiguous base resets the run. The order is `hash64(kmer[z], mask)`, Thomas Wang's 64-bit integer hash with every step masked to 2k bits (invertible, so the hash value is the index key and there are no key collisions); the stored value is `hash<<8 | span`. Ties: when several positions in a window hold the minimal value (necessarily identical canonical k-mers), all of them are emitted ("robust winnowing").
* Presets (`mm_set_opt`, options.c): default / map-ont / lr: k15 w10; map-hifi / lr:hq / map-ccs: k19 w19 (mid_occ clamped to [50, 500]); asm5 / asm10: k19 w19; asm20: k19 w10 (all asm presets clamp mid_occ to [50, 500]); sr: k21 w11 (mid_occ fixed 1000, max_occ 5000); lr:hqae: k25 w51; splice: k15 w5; map-pb: k19 w10 + HPC (`-H`; not studied).
* `-f` (index.c `mm_idx_cal_max_occ` + `mm_mapopt_update`): mid_occ = the (1−f) quantile (f = 2e-4) of per-distinct-minimizer occurrence + 1, clamped to [min_mid_occ, max_mid_occ]; seeds with occurrence > mid_occ are dropped except that `mm_seed_select` keeps up to (streak length/occ_dist) lowest-occurrence seeds (occurrence ≤ 4095) inside high-occurrence streaks (occ_dist = 500; not modelled here).

### Schemes (L1_schemes.py; numba)
| id | scheme | parameters at a preset (k, w) | strand handling |
|---|---|---|---|
| mm2 | minimap2 random minimizer | hash64, robust ties | exact port of `mm_sketch` |
| mm2_1pick | same order, one k-mer per window | — | min over canonical hash |
| rand_mix | random minimizer, strong order | splitmix64 finaliser >>8 | canonical code |
| miniception | closed syncmers first, then random | r = max(4, k−w) | canonical s-mers (symmetric by construction) |
| ddecyc | double decycling (Mykkeltveit classes 0/1/2), random tie-break | as in the reference Rust code | class = min over the two strands |
| greedy | GreedyMini binary order on base>>1 bits, then lex | orders k15-w10 (map-ont, asm20), k15-w19 (map-hifi, k′=15 prefix extension), k16-w11 (sr, k′=16) from OrensteinLab/Minimizers | rank = min over the two strands |
| mod | mod-minimizer | r = 4, t = r + ((k−r) mod w) | canonical t-mer order (see below) |
| modoc | open-closed mod-minimizer | as mod; t-mer order = (open, closed, other) w.r.t. canonical s-mers, s = smallest s ≥ 4 with t−s even | canonical |
| oc | open-closed minimizer | s = smallest s ≥ 5 with k−s even | canonical |
| osync / csync | open / closed syncmers | s chosen to bracket the baseline density (k−s even for osync) | canonical s-mers |

All non-baseline orders use the same splitmix-derived 56-bit random order and the same window logic as `mm_sketch` (including robust ties), except mod-sampling, which picks one k-mer per window.

**Strand.** minimap2 needs a strand-symmetric sketch: the sampled positions of a sequence and of its reverse complement must coincide. Two observations make cheap canonical versions possible:
1. For mod-sampling with canonical t-mer ranks, k − t = w·⌊(k−r)/w⌋ is a multiple of w. A minimal t-mer at window offset x on the forward strand sits at offset w−1+(k−t)−x on the reverse strand, and so selects forward k-mer w−1−((w−1−x) mod w) = x mod w, the same k-mer. Only ties between identical t-mers need a symmetric rule; ties are broken by the hash of the candidate k-mer, and measured symmetry is ≥ 0.998.
2. Open syncmers are symmetric only when the open offset ⌊(k−s)/2⌋ is its own mirror, i.e. k−s is even. With k−s odd, the inclusive canonical rule doubles the open set and destroys the density gain (validation table).

The alternative is a forward-only index with both read strands queried. It keeps the forward density, but doubles query sketching and lookups, and adds the spurious hits of the wrong-orientation sketch. Its index density is reported for comparison.

### Data (public, 1000 Genomes–linked)
* CHM13v2.0 (HPRC bucket `T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz`, md5 06caded0055ec647a69acb708e13beff); chr20 = 66,210,255 bp, no N ([Nurk et al. 2022](https://doi.org/10.1126/science.abj6987)).
* GRCh38 chr20 (UCSC `hg38/chromosomes/chr20.fa.gz`, md5 1d42fbb5014b165373aca053f4f24061), 64,444,167 bp.
* NA19240 HPRC release-2 trio-phased hifiasm 0.19.9 haplotypes (`assemblies_release2_v1.0.index.csv`; [Liao et al. 2023](https://doi.org/10.1038/s41586-023-05896-x)). chr20 was identified from the NCBI assembly reports: paternal `NA19240#1#CM099588.1` (GCA_018503265.2, 66,616,132 bp) and maternal `NA19240#2#CM099626.1` (GCA_018503275.2, 67,338,552 bp). Both were fetched by remote `samtools faidx` from `s3://human-pangenomics/working/HPRC_PLUS/NA19240/assemblies/release2/`.
* Liftover truth: `minimap2 -x asm5 -c --secondary=no` of each haplotype onto CHM13 chr20 (2.31; 414 and 339 alignment records). 98.9 % of paternal bases lift; lifted bases are 99.8 % identical.

### Reads (L1_simulate.py; per-base truth)
* Illumina 2×150: insert N(400, 60); substitutions rising quadratically from 0.1 % to 0.5 % along the read (realised 0.23 %); indels 1e-5. 20,000 pairs from CHM13 and 10,000 pairs from each NA19240 haplotype; mates scored separately.
* HiFi Q30 (realised 0.100 % errors). Composition from [Wenger et al. 2019](https://doi.org/10.1038/s41587-019-0217-9): a mismatch every 13,048 bp, a non-homopolymer indel every 9,669 bp, a homopolymer indel every 477 bp, i.e. 92.0 % homopolymer indels. Homopolymer length changes occur with probability c·(L−1)^1.5 per run, c calibrated on chr20; 65 % are contractions. Per-read rate multiplier lognormal(σ = 0.5); lengths lognormal with mean 16 kb, sd 4 kb. 1,500 reads from CHM13 and 750 per haplotype.
* ONT R10.4.1 Q20 (realised 0.98 %): homopolymer 0.40 % (70 % contractions), substitutions 0.35 %, insertions 0.10 %, deletions 0.15 % per base. This split is illustrative; that homopolymer indels are the residual R10.4 error is only qualitative ([Sereika et al. 2022](https://doi.org/10.1038/s41592-022-01539-7)). Lengths lognormal with mean 30 kb, sd 25 kb. 800 reads from CHM13 and 400 per haplotype.
* ONT, pbsim3 (v3.0.5): ERRHMM-ONT-HQ errhmm, `--accuracy-mean 0.99` ([Ono et al. 2022](https://doi.org/10.1093/nargab/lqac092)). The realised error was 2.95 % (0.84 % substitutions, 0.88 % insertions, 1.24 % deletions), i.e. harder than Q20. Depth 0.36× from CHM13 and 0.18× per haplotype; truth from pbsim3's MAF.
* Assembly scenario: error-free 20-kb tiles of each NA19240 haplotype every 40 kb (1,665 + 1,683 tiles) for map-hifi (= asm5 k, w), asm20 and lr:hqae.

### Metrics (L1_eval.py)
* **Index level** (CHM13 chr20; GRCh38 chr20 for comparison): entries; distinct keys; density (entries/valid k-mers); density factor d·(w+1); mid_occ by minimap2's rule; fraction of entries above mid_occ ("high-frequency"); largest gap between consecutive sampled positions.
* **Per read.** Each read seed is looked up in the index. Its truth is the lifted position of the two end bases of the k-mer: identity for CHM13 reads, the PAF liftover for NA19240 reads. A seed is a true-origin match if any index hit has the expected relative strand and lies within 20 bp of the expected end (k+10 if one end is unliftable).
  * Seeds above mid_occ are dropped, as with minimap2's hard cutoff.
  * Reported per read: true matches, fraction of reads with ≥ 1 and ≥ 3 true matches, anchors (Σ occurrences of retained seeds), and the longest read stretch not covered by a true match.
  * Reads whose truth is undefined (< 50 % of baseline seeds liftable) are excluded for all schemes.
  * Strata, defined on the baseline: unique reads (< 50 % of seeds above mid_occ) and repeat reads (the rest).
* **Statistics.** Per-read metrics use a paired bootstrap over reads (B = 1,000; the same resampled reads for scheme and baseline), giving 95 % CIs for % change of means and for differences in fractions. Index entries use a bootstrap over 66 one-megabase reference blocks.
* **Equal-density curves.** w is varied for windowed schemes and s for syncmers, at the preset k, on NA19240 reads.

## Results

### Validation
The numba port of `mm_sketch` reproduces minimap2 2.31's own `sketch.c` bit-for-bit, compiled into a 30-line C harness, on all 66.2 Mb of CHM13 chr20. It matches at sr, map-ont, map-hifi, asm20, lr:hqae and splice, i.e. 2.57–22.6 M (x, y) pairs per preset. The patched C routine (`--mod-oc`) and the Python `modoc` select identical position sets at all five presets (2.23–9.39 M positions).

Published densities reproduce on a 10-Mb i.i.d. random string:

| Scheme (w, k, s) | Measured | Exact / theory |
|---|---|---|
| Miniception (5, 11, 6) | 0.2934 | 0.2929 |
| Open-closed minimizer (5, 11, 6) | 0.2868 | 0.2864 |
| Mod-minimizer (10, 15) | 0.1446 | 0.1429 (assumes no duplicate t-mers) |

Strand symmetry, measured as the fraction of positions shared by sketch(S) and mirror(sketch(rc S)), is 1.0000 for every order- and syncmer-based scheme and 0.9976–0.9990 for mod-sampling, where tied identical t-mers break the symmetry.

### Finding 1: minimap2's hash is a weak random order at k ≤ 15
`hash64(key, mask)` masks every step to 2k bits. For 2k ≤ 31 the `key<<31` step is a no-op and the `>>24`/`>>28` mixes act on few or no bits. On random sequence, minimap2-exact density at map-ont (w10, k15) is 0.18655. A strong order gives 0.18196, and a width-adapted invertible 2k-bit mixer (`hash_w`, bijective by construction) gives 0.18180. The ideal 2/(w+1) is 0.18182. At splice (w5, k15) the values are 0.34007 vs 0.33342. At k = 19, 21 and 25 there is no difference.

On CHM13 chr20 the strong order alone cuts map-ont index entries by 2.8 % [2.55, 3.05] and true seeds per read by 2.9 %, i.e. in proportion to density. It is a free, format-compatible fix whenever the index is rebuilt.

### Finding 2: canonicalization decides which low-density schemes survive
minimap2 needs strand-symmetric sketches. Taking min(rank(x), rank(rc x)) erodes most of the gain of orders optimised on one strand. On random sequence, GreedyMini at w10/k15 goes from 0.1344 (forward) to 0.1671 (canonical), and double decycling at w19/k19 from 0.0790 to 0.0870. Open-closed schemes are fully symmetric at no cost once s is chosen with k−s (resp. t−s) even: w19/k19, s = 5 gives 0.0802 forward vs 0.0808 canonical. With k−s odd it gives 0.0799 vs 0.0988. Mod-sampling is symmetric by construction, because k−t is a multiple of w.

### Finding 3: the open-closed mod-minimizer gives the largest index reduction at every preset
The open-closed mod-minimizer cuts CHM13 chr20 index entries at every preset (Table 1):

| Preset | modoc | mod-minimizer | Next best |
|---|---|---|---|
| sr | −26.1 % | −21.8 % | — |
| map-ont | −24.6 % | −22.3 % | OC −20.9 % |
| map-hifi / lr:hq / asm5 | −18.8 % | −0.1 % (degenerates, t = k) | OC −18.4 % |
| asm20 | −24.6 % | −21.1 % | — |
| lr:hqae | −13.1 % | −0.8 % | — |

Density factors fall from 2.00–2.07 to 1.48–1.75; the simple approximation ⌈(w+k)/w⌉/(w+k) of the forward-scheme lower bound gives 1.05–1.37 (the near-tight bound is higher unless k ≡ 1 mod w; dashed line in panel a). GRCh38 chr20 reproduces each reduction to within 0.2 pp. Miniception, double decycling and GreedyMini (with naive canonicalization) reach −5 % to −16 %. The high-frequency share of entries changes little: at map-hifi it is 8.6 % for minimap2, 9.0 % for modoc and 9.4 % for OC.

### Finding 4: in unique sequence the density cut costs no read sensitivity
True-origin seeds per read fall almost exactly in proportion to density. For NA19240 ONT reads modoc gives −23.5 % [−23.6, −23.4], against −24.6 % index entries. Anchors per read fall by 22.6 % (sr), 22.1 % (ONT model), 22.5 % (pbsim3 ONT), 13.8 % (HiFi), 14.3 % (asm5 tiles), 21.0 % (asm20 tiles) and 8.5 % (lr:hqae tiles).

The fraction of reads with ≥ 3 true seeds is unchanged for unique-sequence long reads and tiles. Unique reads are 94.2–99.4 % of reads, stratum defined on the baseline as < 50 % of seeds above mid_occ. All schemes give 100.00 %, Δ 0.00 pp. For CHM13-derived reads, where only sequencing error applies, every long-read scheme gives 100 %.

For unique Illumina reads, ≥ 1 true seed changes by −0.01 pp [−0.02, +0.01] and ≥ 3 by −0.11 pp [−0.15, −0.07]. The 95th-percentile longest seedless gap improves from 35 to 32 bp, because mod-sampling spreads samples evenly. For unique long reads the p95 gap changes by +3 % (ONT, 318 vs 309 bp), −2 % (HiFi) and +10 % (asm5 tiles, 174 vs 158 bp).

### Finding 5: the risk is concentrated in repeat-dominated reads
Repeat-dominated reads are 0.6–5.8 % of NA19240 reads and tiles, mostly satellites. There, the outcome depends on a handful of seeds below mid_occ, and any change of order perturbs it. The strong-hash random minimizer alone moves ≥ 1 true seed by −10 pp (HiFi) and true seeds by +31 % (lr:hqae tiles). The low-density schemes lose more; for modoc, ≥ 1 true seed falls by:

| Preset | Repeat reads (n) | Δ ≥ 1 true seed, pp [95 % CI] |
|---|---|---|
| map-ont | 40 | −20 [−32.5, −7.5] |
| map-hifi | 80 | −16 [−25, −9] |
| asm5 tiles | 191 | −12.6 |
| lr:hqae tiles | 189 | −4.8 |

This metric uses minimap2's hard mid_occ cut-off and does not model `mm_seed_select`, which rescues up to one low-occurrence seed per 500 bp of high-occurrence streak. In the real-aligner test below, mapping accuracy was non-inferior.

### Finding 6: equal density favours context-free sampling for short reads only
At equal index density, closed syncmers give Illumina reads +9.7 % [9.5, 9.9] true seeds per read (+0.8 % entries, s = 10) and a shorter p95 gap (31 vs 35 bp). Part of this is the read-end effect: minimizers need a complete window, while syncmers sample every qualifying k-mer up to the read ends. For long reads, open/closed syncmers add true seeds roughly in proportion to the entries they add:

| Preset | Scheme | Entries | True seeds | Anchors |
|---|---|---|---|---|
| map-ont | open syncmer, s = 11 | +6.6 % | +8.4 % | +15.8 % |
| map-hifi | open syncmer, s = 11 | +11.5 % | +11.7 % | +22.0 % |

Anchors grow about twice as fast because context-free sampling also samples every copy in tandem repeats. Syncmers therefore help short-read sensitivity but are not an anchor-efficient choice for long reads. 

Equal-density curves (panels d–f; [L1_equal_density.tsv](L1_equal_density.tsv), [L1_curves.tsv](L1_curves.tsv)) make this explicit. Interpolated at minimap2's preset density, Illumina reads gain +7.8 % (mod-minimizer), +8.0 % (open-closed mod-minimizer), +9.3 % (open syncmer) and +8.9 % (closed syncmer) true seeds per kb, with p95 gaps 16 %, 16 %, 15 % and 7 % shorter. ONT reads gain only +1.0 to +1.5 %, while anchors per kb rise 6.8–8.9 %, and HiFi reads +0.0 to +0.5 %. For long reads all schemes are equally sensitive per indexed k-mer, so their value is fewer entries and anchors, not more true seeds at equal index size.

### Finding 7: the patched minimap2 is smaller and faster at equal accuracy
`minimap2_sampling.patch` adds `--mod-oc`: index flag 0x8 plus a new .mmi magic `MMI\3`. Old binaries therefore do not recognise the file as an index. In a test, stock 2.31 parsed it as FASTA and produced no alignments (exit 0, 0/400 reads mapped) rather than confident wrong ones. This is a no-output failure, not an explicit error message. Without the flag, the patched binary writes indexes byte-identical to stock. On the CHM13 chr20 index with 1 thread (Table 4):

| Preset | .mmi size | Mapping CPU time | Wrong mappings, MAPQ ≥ 1 | Index build CPU |
|---|---|---|---|---|
| sr | −21.9 % | −20 % (5.23 → 4.18 s) | 0.319 % → 0.271 % | 1.45× |
| map-ont | −21.0 % | −30 % (37.7 → 26.3 s) | 0.259 % → 0.259 % | 1.30× |
| map-hifi | −14.3 % | −12 % (14.9 → 13.1 s) | 0.345 % → 0.138 % | 1.85× |

MAPQ-60 counts are unchanged and all reads map with both binaries. Index construction is slower because the routine is unoptimised scalar code. 

**Whole genome (CHM13v2.0, 1 thread, split index; [L1_wholegenome.tsv](L1_wholegenome.tsv)).** At map-hifi (`-I 1G`) the `--mod-oc` index is 14.2 % smaller (5,141,230,810 vs 5,993,367,034 bytes) and its peak build RSS is 33 % lower (5.19 vs 7.70 GB). Construction takes 1.9× the CPU (151.3 vs 78.8 s). Genome-wide mapping of the 2,942 chr20 HiFi reads (`--split-prefix`) gave 2,880 vs 2,847 reads at MAPQ 60 (1 vs 0 wrong), and wrong mappings at MAPQ ≥ 1 of 0.310 % vs 0.482 %. Mapping CPU was equal (32.1 vs 31.4 s) because loading the four-part index dominates for so few reads. At map-ont (`-I 500M`) the index is 22.1 % smaller (6,877,839,462 vs 8,825,362,630 bytes), with 1.31× build CPU and a higher peak build RSS (8.71 vs 7.18 GB).

## Tables

**Table 1. Index statistics, CHM13v2.0 chr20.** Δentries 95 % CI from a bootstrap over 66 one-megabase blocks; the GRCh38 column gives the same comparison on GRCh38 chr20; high-freq = entries whose key occurs more than mid_occ times.

| preset (k,w) | scheme | s | density | d·(w+1) | entries | Δentries vs mm2 (95% CI) | mid_occ | high-freq entries | GRCh38 Δentries |
|---|---|---|---|---|---|---|---|---|---|
| sr (21,11) | minimap2 2.31 (baseline) |  | 0.1672 | 2.006 | 11,069,512 | +0.0% | 1000 | 2.7% | +0.0% |
| sr (21,11) | same order, 1 k-mer/window |  | 0.1670 | 2.004 | 11,058,486 | -0.1% [-0.11, -0.09] | 1000 | 2.6% | -0.1% |
| sr (21,11) | random, strong hash |  | 0.1672 | 2.006 | 11,069,954 | +0.0% [-0.12, 0.10] | 1000 | 2.7% | -0.0% |
| sr (21,11) | miniception |  | 0.1429 | 1.714 | 9,458,460 | -14.6% [-14.67, -14.47] | 1000 | 3.0% | -14.5% |
| sr (21,11) | double decycling |  | 0.1486 | 1.783 | 9,839,851 | -11.1% [-11.18, -11.01] | 1000 | 2.9% | -11.2% |
| sr (21,11) | GreedyMini (min over strands) |  | 0.1427 | 1.712 | 9,447,427 | -14.7% [-14.91, -14.45] | 1000 | 2.9% | -14.6% |
| sr (21,11) | mod-minimizer |  | 0.1307 | 1.569 | 8,656,676 | -21.8% [-21.91, -21.71] | 1000 | 2.9% | -21.8% |
| sr (21,11) | open-closed mod-minimizer |  | 0.1235 | 1.482 | 8,175,951 | -26.1% [-26.18, -26.09] | 1000 | 3.1% | -26.2% |
| sr (21,11) | open-closed minimizer |  | 0.1346 | 1.615 | 8,912,424 | -19.5% [-19.61, -19.37] | 1000 | 3.3% | -19.5% |
| sr (21,11) | open syncmer | 15 | 0.1436 | 1.723 | 9,507,435 | -14.1% [-14.22, -13.98] | 1000 | 3.5% | -14.2% |
| sr (21,11) | closed syncmer | 10 | 0.1685 | 2.022 | 11,159,151 | +0.8% [0.51, 1.19] | 1000 | 3.5% | +0.7% |
| map-ont (15,10) | minimap2 2.31 (baseline) |  | 0.1881 | 2.069 | 12,452,943 | +0.0% | 201 | 10.6% | +0.0% |
| map-ont (15,10) | same order, 1 k-mer/window |  | 0.1877 | 2.065 | 12,426,845 | -0.2% [-0.25, -0.18] | 201 | 10.4% | -0.2% |
| map-ont (15,10) | minimap2 order, width-adapted invertible hash |  | 0.1828 | 2.011 | 12,106,122 | -2.8% [-2.93, -2.59] | 196 | 10.6% |  |
| map-ont (15,10) | random, strong hash |  | 0.1828 | 2.010 | 12,101,017 | -2.8% [-3.05, -2.55] | 200 | 10.6% | -2.9% |
| map-ont (15,10) | miniception |  | 0.1581 | 1.739 | 10,469,227 | -15.9% [-16.05, -15.76] | 213 | 11.2% | -16.0% |
| map-ont (15,10) | double decycling |  | 0.1589 | 1.747 | 10,517,721 | -15.5% [-15.60, -15.47] | 205 | 11.0% | -15.5% |
| map-ont (15,10) | GreedyMini (min over strands) |  | 0.1671 | 1.838 | 11,064,696 | -11.1% [-11.22, -11.07] | 206 | 10.6% | -11.2% |
| map-ont (15,10) | mod-minimizer |  | 0.1461 | 1.607 | 9,674,088 | -22.3% [-22.50, -22.08] | 200 | 10.7% | -22.4% |
| map-ont (15,10) | open-closed mod-minimizer |  | 0.1418 | 1.560 | 9,388,774 | -24.6% [-24.72, -24.50] | 198 | 10.5% | -24.6% |
| map-ont (15,10) | open-closed minimizer |  | 0.1488 | 1.636 | 9,850,006 | -20.9% [-21.02, -20.76] | 208 | 11.2% | -21.0% |
| map-ont (15,10) | open syncmer | 11 | 0.2005 | 2.206 | 13,277,158 | +6.6% [6.46, 6.81] | 211 | 11.0% | +6.6% |
| map-ont (15,10) | closed syncmer | 5 | 0.1904 | 2.094 | 12,605,670 | +1.2% [1.10, 1.37] | 227 | 11.6% | +1.2% |
| map-hifi (19,19) | minimap2 2.31 (baseline) |  | 0.1005 | 2.010 | 6,652,482 | +0.0% | 146 | 8.6% | +0.0% |
| map-hifi (19,19) | same order, 1 k-mer/window |  | 0.1003 | 2.005 | 6,638,460 | -0.2% [-0.24, -0.19] | 146 | 8.4% | -0.2% |
| map-hifi (19,19) | random, strong hash |  | 0.1006 | 2.013 | 6,662,397 | +0.1% [-0.02, 0.45] | 146 | 8.6% | +0.0% |
| map-hifi (19,19) | miniception |  | 0.0912 | 1.824 | 6,037,405 | -9.2% [-9.55, -8.92] | 149 | 8.9% | -9.3% |
| map-hifi (19,19) | double decycling |  | 0.0876 | 1.751 | 5,797,392 | -12.9% [-13.02, -12.61] | 148 | 9.0% | -13.0% |
| map-hifi (19,19) | GreedyMini (min over strands) |  | 0.0887 | 1.773 | 5,870,367 | -11.8% [-12.01, -11.44] | 156 | 9.1% | -11.9% |
| map-hifi (19,19) | mod-minimizer |  | 0.1004 | 2.007 | 6,645,850 | -0.1% [-0.20, 0.02] | 145 | 8.4% | -0.1% |
| map-hifi (19,19) | open-closed mod-minimizer |  | 0.0816 | 1.632 | 5,402,121 | -18.8% [-18.92, -18.61] | 153 | 9.0% | -18.9% |
| map-hifi (19,19) | open-closed minimizer |  | 0.0820 | 1.640 | 5,427,820 | -18.4% [-18.63, -18.04] | 154 | 9.4% | -18.6% |
| map-hifi (19,19) | open syncmer | 11 | 0.1121 | 2.241 | 7,419,462 | +11.5% [11.09, 11.97] | 156 | 9.4% | +11.4% |
| asm20 (19,10) | minimap2 2.31 (baseline) |  | 0.1820 | 2.002 | 12,050,520 | +0.0% | 146 | 8.6% | +0.0% |
| asm20 (19,10) | same order, 1 k-mer/window |  | 0.1818 | 2.000 | 12,035,444 | -0.1% [-0.14, -0.11] | 146 | 8.4% | -0.1% |
| asm20 (19,10) | random, strong hash |  | 0.1827 | 2.009 | 12,093,272 | +0.4% [0.11, 0.65] | 149 | 8.8% | +0.2% |
| asm20 (19,10) | miniception |  | 0.1571 | 1.728 | 10,401,977 | -13.7% [-13.99, -13.31] | 156 | 9.3% | -13.8% |
| asm20 (19,10) | double decycling |  | 0.1618 | 1.780 | 10,712,684 | -11.1% [-11.26, -10.92] | 147 | 8.8% | -11.1% |
| asm20 (19,10) | GreedyMini (min over strands) |  | 0.1555 | 1.710 | 10,293,589 | -14.6% [-14.63, -14.53] | 151 | 8.9% | -14.6% |
| asm20 (19,10) | mod-minimizer |  | 0.1436 | 1.580 | 9,509,469 | -21.1% [-21.30, -20.80] | 154 | 8.9% | -21.2% |
| asm20 (19,10) | open-closed mod-minimizer |  | 0.1372 | 1.509 | 9,085,652 | -24.6% [-24.75, -24.41] | 151 | 8.9% | -24.7% |
| asm20 (19,10) | open-closed minimizer |  | 0.1484 | 1.632 | 9,822,423 | -18.5% [-18.59, -18.37] | 150 | 9.2% | -18.6% |
| asm20 (19,10) | open syncmer | 15 | 0.2001 | 2.201 | 13,247,778 | +9.9% [9.84, 10.05] | 150 | 9.0% | +9.9% |
| asm20 (19,10) | closed syncmer | 9 | 0.1844 | 2.028 | 12,206,646 | +1.3% [0.98, 1.66] | 162 | 9.6% | +1.2% |
| lr:hqae (25,51) | minimap2 2.31 (baseline) |  | 0.0388 | 2.018 | 2,569,321 | +0.0% | 109 | 6.7% | +0.0% |
| lr:hqae (25,51) | same order, 1 k-mer/window |  | 0.0387 | 2.015 | 2,565,575 | -0.1% [-0.24, -0.09] | 109 | 6.6% | -0.1% |
| lr:hqae (25,51) | random, strong hash |  | 0.0386 | 2.007 | 2,555,406 | -0.5% [-0.90, -0.26] | 110 | 6.3% | -0.5% |
| lr:hqae (25,51) | miniception |  | 0.0369 | 1.920 | 2,444,741 | -4.8% [-5.14, -4.60] | 114 | 6.8% | -4.9% |
| lr:hqae (25,51) | double decycling |  | 0.0359 | 1.865 | 2,374,229 | -7.6% [-8.08, -7.28] | 111 | 7.0% | -7.4% |
| lr:hqae (25,51) | mod-minimizer |  | 0.0385 | 2.002 | 2,549,287 | -0.8% [-1.27, -0.41] | 110 | 6.1% | -0.6% |
| lr:hqae (25,51) | open-closed mod-minimizer |  | 0.0337 | 1.754 | 2,233,490 | -13.1% [-13.27, -12.90] | 111 | 7.3% | -13.0% |
| lr:hqae (25,51) | open-closed minimizer |  | 0.0339 | 1.761 | 2,241,906 | -12.7% [-13.03, -12.37] | 112 | 7.6% | -12.8% |

**Table 2. Seed metrics for NA19240-derived reads and tiles** (sequencing error plus real divergence; hap = NA19240 pat+mat reads, pbhap = pbsim3 ONT-HQ reads, asm = error-free 20-kb tiles). Δ values are vs minimap2 with paired-bootstrap 95 % CIs over reads (B = 1,000).

| preset | reads | scheme | true seeds/read | Δ true (95% CI) | Δ anchors (95% CI) | ≥1 true | ≥3 true | Δ≥3 pp (95% CI) | longest gap mean / p95 (bp) |
|---|---|---|---|---|---|---|---|---|---|
| sr | hap (n=39338) | minimap2 2.31 (baseline) | 17.7 | — | — | 99.97% | 99.83% | — | 16 / 35 |
| sr | hap (n=39338) | random, strong hash | 17.7 | +0.1% [-0.10, 0.23] | -0.5% [-1.35, 0.41] | 99.96% | 99.80% | -0.04 [-0.081, 0.003] | 16 / 36 |
| sr | hap (n=39338) | miniception | 15.1 | -14.3% [-14.44, -14.17] | -13.9% [-14.66, -13.15] | 99.96% | 99.77% | -0.07 [-0.112, -0.025] | 17 / 37 |
| sr | hap (n=39338) | double decycling | 15.7 | -11.0% [-11.11, -10.81] | -13.0% [-13.66, -12.20] | 99.98% | 99.77% | -0.06 [-0.104, -0.020] | 17 / 36 |
| sr | hap (n=39338) | GreedyMini (min over strands) | 15.1 | -14.4% [-14.50, -14.23] | -15.0% [-15.69, -14.25] | 99.96% | 99.72% | -0.11 [-0.160, -0.069] | 17 / 36 |
| sr | hap (n=39338) | mod-minimizer | 14.7 | -16.7% [-16.86, -16.60] | -20.3% [-21.09, -19.50] | 99.95% | 99.68% | -0.15 [-0.203, -0.099] | 12 / 33 |
| sr | hap (n=39338) | open-closed mod-minimizer | 13.9 | -21.1% [-21.25, -21.02] | -22.6% [-23.27, -21.91] | 99.95% | 99.64% | -0.19 [-0.247, -0.142] | 12 / 33 |
| sr | hap (n=39338) | open-closed minimizer | 14.2 | -19.4% [-19.48, -19.24] | -24.6% [-25.27, -23.98] | 99.97% | 99.70% | -0.13 [-0.178, -0.086] | 17 / 36 |
| sr | hap (n=39338) | open syncmer | 16.6 | -6.3% [-6.41, -6.12] | -9.3% [-10.04, -8.50] | 99.98% | 99.80% | -0.04 [-0.074, 0.005] | 10 / 31 |
| sr | hap (n=39338) | closed syncmer | 19.4 | +9.7% [9.55, 9.91] | +10.4% [9.43, 11.36] | 99.97% | 99.85% | +0.02 [-0.020, 0.056] | 11 / 31 |
| map-ont | hap (n=795) | minimap2 2.31 (baseline) | 4060.5 | — | — | 99.75% | 99.75% | — | 568 / 2422 |
| map-ont | hap (n=795) | minimap2 order, width-adapted invertible hash | 3947.7 | -2.8% [-2.85, -2.70] | -4.9% [-5.56, -4.29] | 99.75% | 99.62% | -0.13 [-0.377, 0.000] | 567 / 2438 |
| map-ont | hap (n=795) | random, strong hash | 3944.6 | -2.9% [-2.93, -2.78] | -3.8% [-4.24, -3.29] | 99.75% | 99.37% | -0.38 [-0.881, 0.000] | 551 / 2409 |
| map-ont | hap (n=795) | miniception | 3437.6 | -15.3% [-15.42, -15.26] | -8.8% [-9.66, -8.10] | 99.75% | 98.87% | -0.88 [-1.635, -0.252] | 570 / 2336 |
| map-ont | hap (n=795) | double decycling | 3451.3 | -15.0% [-15.08, -14.92] | -13.0% [-13.59, -12.41] | 99.75% | 98.99% | -0.75 [-1.384, -0.252] | 564 / 2348 |
| map-ont | hap (n=795) | GreedyMini (min over strands) | 3654.6 | -10.0% [-10.09, -9.90] | -6.1% [-6.91, -5.45] | 99.75% | 98.99% | -0.75 [-1.384, -0.252] | 583 / 2261 |
| map-ont | hap (n=795) | mod-minimizer | 3188.3 | -21.5% [-21.56, -21.41] | -18.3% [-18.93, -17.78] | 99.75% | 98.99% | -0.75 [-1.384, -0.252] | 586 / 2482 |
| map-ont | hap (n=795) | open-closed mod-minimizer | 3107.2 | -23.5% [-23.57, -23.39] | -22.1% [-22.89, -21.40] | 98.74% | 97.61% | -2.14 [-3.145, -1.258] | 727 / 2499 |
| map-ont | hap (n=795) | open-closed minimizer | 3233.2 | -20.4% [-20.44, -20.30] | -15.7% [-16.23, -15.24] | 99.75% | 98.99% | -0.75 [-1.384, -0.252] | 586 / 2416 |
| map-ont | hap (n=795) | open syncmer | 4400.3 | +8.4% [8.25, 8.49] | +15.8% [15.12, 16.40] | 99.75% | 99.62% | -0.13 [-0.377, 0.000] | 588 / 2391 |
| map-ont | hap (n=795) | closed syncmer | 4138.9 | +1.9% [1.82, 2.05] | +21.4% [20.43, 22.50] | 99.75% | 99.12% | -0.63 [-1.258, -0.126] | 553 / 2307 |
| map-ont | pbhap (n=779) | minimap2 2.31 (baseline) | 3884.2 | — | — | 99.36% | 98.97% | — | 756 / 3889 |
| map-ont | pbhap (n=779) | minimap2 order, width-adapted invertible hash | 3776.0 | -2.8% [-2.87, -2.70] | -5.6% [-6.59, -4.75] | 99.23% | 98.84% | -0.13 [-0.517, 0.257] | 796 / 4129 |
| map-ont | pbhap (n=779) | random, strong hash | 3774.5 | -2.8% [-2.90, -2.74] | -3.6% [-4.21, -2.99] | 99.36% | 98.59% | -0.39 [-0.902, 0.128] | 807 / 3886 |
| map-ont | pbhap (n=779) | miniception | 3293.2 | -15.2% [-15.34, -15.10] | -8.7% [-9.31, -8.03] | 99.36% | 98.46% | -0.51 [-1.027, -0.128] | 735 / 3738 |
| map-ont | pbhap (n=779) | double decycling | 3300.8 | -15.0% [-15.11, -14.92] | -12.6% [-13.49, -11.73] | 99.49% | 98.72% | -0.26 [-0.770, 0.257] | 773 / 4105 |
| map-ont | pbhap (n=779) | GreedyMini (min over strands) | 3492.3 | -10.1% [-10.20, -9.98] | -6.5% [-7.02, -5.98] | 99.10% | 98.72% | -0.26 [-0.770, 0.257] | 801 / 4030 |
| map-ont | pbhap (n=779) | mod-minimizer | 3050.2 | -21.5% [-21.58, -21.37] | -18.7% [-19.26, -18.10] | 99.23% | 98.72% | -0.26 [-0.899, 0.385] | 819 / 3886 |
| map-ont | pbhap (n=779) | open-closed mod-minimizer | 2971.8 | -23.5% [-23.60, -23.39] | -22.5% [-23.45, -21.65] | 98.20% | 97.95% | -1.03 [-1.797, -0.257] | 1170 / 5279 |
| map-ont | pbhap (n=779) | open-closed minimizer | 3093.4 | -20.4% [-20.46, -20.25] | -15.9% [-16.38, -15.37] | 99.23% | 98.59% | -0.39 [-1.027, 0.128] | 823 / 4028 |
| map-ont | pbhap (n=779) | open syncmer | 4212.3 | +8.4% [8.26, 8.66] | +16.0% [15.04, 16.87] | 99.23% | 99.10% | +0.13 [-0.257, 0.513] | 824 / 4030 |
| map-ont | pbhap (n=779) | closed syncmer | 3968.9 | +2.2% [1.99, 2.37] | +22.1% [21.21, 22.93] | 99.23% | 98.84% | -0.13 [-0.770, 0.385] | 744 / 2951 |
| map-hifi | hap (n=1482) | minimap2 2.31 (baseline) | 1376.2 | — | — | 99.46% | 98.99% | — | 454 / 2077 |
| map-hifi | hap (n=1482) | random, strong hash | 1377.7 | +0.1% [0.01, 0.22] | +0.3% [-0.46, 1.07] | 98.92% | 98.52% | -0.47 [-0.877, -0.135] | 518 / 2222 |
| map-hifi | hap (n=1482) | miniception | 1250.2 | -9.2% [-9.24, -9.06] | -6.3% [-7.05, -5.60] | 98.65% | 98.18% | -0.81 [-1.282, -0.405] | 598 / 2426 |
| map-hifi | hap (n=1482) | double decycling | 1200.1 | -12.8% [-12.87, -12.72] | -8.6% [-9.45, -7.71] | 98.92% | 98.18% | -0.81 [-1.282, -0.405] | 536 / 2211 |
| map-hifi | hap (n=1482) | GreedyMini (min over strands) | 1213.9 | -11.8% [-11.88, -11.71] | -5.7% [-7.07, -4.32] | 98.58% | 97.44% | -1.55 [-2.227, -0.945] | 607 / 2420 |
| map-hifi | hap (n=1482) | mod-minimizer | 1378.3 | +0.2% [0.05, 0.26] | -0.3% [-1.08, 0.48] | 98.92% | 98.52% | -0.47 [-0.810, -0.135] | 519 / 2222 |
| map-hifi | hap (n=1482) | open-closed mod-minimizer | 1120.5 | -18.6% [-18.65, -18.50] | -13.8% [-14.35, -13.14] | 98.58% | 97.37% | -1.62 [-2.362, -1.012] | 600 / 2431 |
| map-hifi | hap (n=1482) | open-closed minimizer | 1119.9 | -18.6% [-18.69, -18.55] | -13.2% [-13.90, -12.52] | 98.58% | 97.37% | -1.62 [-2.227, -1.012] | 599 / 2431 |
| map-hifi | hap (n=1482) | open syncmer | 1537.9 | +11.7% [11.64, 11.85] | +22.0% [20.95, 23.10] | 98.52% | 97.57% | -1.42 [-2.092, -0.877] | 602 / 2222 |
| map-hifi | asm (n=3315) | minimap2 2.31 (baseline) | 1749.0 | — | — | 99.31% | 98.73% | — | 590 / 2888 |
| map-hifi | asm (n=3315) | random, strong hash | 1750.0 | +0.1% [-0.00, 0.12] | -0.2% [-0.66, 0.32] | 98.94% | 98.28% | -0.45 [-0.724, -0.210] | 651 / 3118 |
| map-hifi | asm (n=3315) | miniception | 1586.2 | -9.3% [-9.36, -9.26] | -6.5% [-6.93, -5.98] | 98.70% | 98.07% | -0.66 [-0.965, -0.362] | 723 / 3207 |
| map-hifi | asm (n=3315) | double decycling | 1521.8 | -13.0% [-13.04, -12.94] | -8.6% [-9.08, -8.18] | 98.85% | 98.07% | -0.66 [-0.935, -0.392] | 667 / 3050 |
| map-hifi | asm (n=3315) | GreedyMini (min over strands) | 1538.7 | -12.0% [-12.08, -11.97] | -6.5% [-7.22, -5.62] | 98.64% | 97.62% | -1.12 [-1.478, -0.754] | 727 / 3078 |
| map-hifi | asm (n=3315) | mod-minimizer | 1750.5 | +0.1% [0.03, 0.15] | -0.8% [-1.29, -0.33] | 98.91% | 98.28% | -0.45 [-0.724, -0.211] | 652 / 3118 |
| map-hifi | asm (n=3315) | open-closed mod-minimizer | 1420.0 | -18.8% [-18.85, -18.77] | -14.3% [-14.75, -13.95] | 98.58% | 97.53% | -1.21 [-1.599, -0.845] | 730 / 3288 |
| map-hifi | asm (n=3315) | open-closed minimizer | 1419.4 | -18.8% [-18.89, -18.81] | -13.9% [-14.31, -13.48] | 98.58% | 97.53% | -1.21 [-1.599, -0.845] | 729 / 3288 |
| map-hifi | asm (n=3315) | open syncmer | 1947.0 | +11.3% [11.26, 11.38] | +20.9% [20.24, 21.53] | 98.58% | 97.71% | -1.03 [-1.418, -0.694] | 727 / 3279 |
| asm20 | asm (n=3315) | minimap2 2.31 (baseline) | 3180.5 | — | — | 98.88% | 98.49% | — | 632 / 2961 |
| asm20 | asm (n=3315) | random, strong hash | 3183.2 | +0.1% [0.04, 0.13] | +1.3% [0.95, 1.57] | 98.79% | 98.49% | +0.00 [-0.181, 0.181] | 646 / 2990 |
| asm20 | asm (n=3315) | miniception | 2727.4 | -14.2% [-14.28, -14.21] | -8.0% [-8.34, -7.70] | 98.73% | 98.40% | -0.09 [-0.302, 0.121] | 652 / 2996 |
| asm20 | asm (n=3315) | double decycling | 2822.2 | -11.3% [-11.30, -11.22] | -11.9% [-12.32, -11.43] | 98.67% | 98.25% | -0.24 [-0.452, -0.060] | 715 / 3121 |
| asm20 | asm (n=3315) | GreedyMini (min over strands) | 2713.3 | -14.7% [-14.73, -14.65] | -12.2% [-12.61, -11.86] | 98.61% | 98.31% | -0.18 [-0.362, 0.000] | 716 / 3065 |
| asm20 | asm (n=3315) | mod-minimizer | 2504.0 | -21.3% [-21.30, -21.24] | -16.5% [-16.89, -16.18] | 98.67% | 98.49% | +0.00 [-0.151, 0.151] | 652 / 3007 |
| asm20 | asm (n=3315) | open-closed mod-minimizer | 2393.9 | -24.7% [-24.76, -24.70] | -21.0% [-21.30, -20.71] | 98.67% | 98.28% | -0.21 [-0.452, 0.000] | 712 / 3069 |
| asm20 | asm (n=3315) | open-closed minimizer | 2580.6 | -18.9% [-18.90, -18.83] | -14.7% [-15.03, -14.39] | 98.61% | 98.25% | -0.24 [-0.452, -0.060] | 713 / 3005 |
| asm20 | asm (n=3315) | open syncmer | 3493.4 | +9.8% [9.79, 9.89] | +14.1% [13.61, 14.55] | 98.61% | 98.40% | -0.09 [-0.271, 0.090] | 712 / 3156 |
| asm20 | asm (n=3315) | closed syncmer | 3196.2 | +0.5% [0.45, 0.54] | +13.6% [13.06, 14.09] | 98.67% | 98.49% | +0.00 [-0.152, 0.181] | 653 / 2991 |
| lr:hqae | asm (n=3315) | minimap2 2.31 (baseline) | 676.1 | — | — | 98.58% | 97.22% | — | 817 / 4756 |
| lr:hqae | asm (n=3315) | random, strong hash | 674.5 | -0.2% [-0.33, -0.13] | +5.2% [4.12, 6.34] | 98.79% | 98.04% | +0.81 [0.452, 1.176] | 733 / 3865 |
| lr:hqae | asm (n=3315) | miniception | 643.9 | -4.8% [-4.86, -4.67] | -1.4% [-2.22, -0.57] | 98.40% | 97.53% | +0.30 [0.000, 0.633] | 803 / 3832 |
| lr:hqae | asm (n=3315) | double decycling | 624.5 | -7.6% [-7.72, -7.54] | -5.1% [-5.77, -4.41] | 98.76% | 97.86% | +0.63 [0.332, 0.936] | 762 / 4594 |
| lr:hqae | asm (n=3315) | mod-minimizer | 675.3 | -0.1% [-0.22, -0.02] | +4.8% [3.64, 5.94] | 98.82% | 98.04% | +0.81 [0.452, 1.176] | 732 / 3865 |
| lr:hqae | asm (n=3315) | open-closed mod-minimizer | 589.0 | -12.9% [-12.96, -12.80] | -8.5% [-9.33, -7.59] | 98.31% | 96.95% | -0.27 [-0.513, -0.030] | 843 / 4879 |
| lr:hqae | asm (n=3315) | open-closed minimizer | 588.3 | -13.0% [-13.07, -12.90] | -8.1% [-9.00, -7.33] | 98.31% | 96.95% | -0.27 [-0.543, 0.000] | 844 / 4879 |

**Table 3. Seed metrics for CHM13-derived reads** (sequencing error only; pbself = pbsim3 ONT-HQ).

| preset | reads | scheme | true seeds/read | Δ true (95% CI) | Δ anchors (95% CI) | ≥1 true | ≥3 true | Δ≥3 pp (95% CI) | longest gap mean / p95 (bp) |
|---|---|---|---|---|---|---|---|---|---|
| sr | self (n=40000) | minimap2 2.31 (baseline) | 18.4 | — | — | 100.00% | 100.00% | — | 14 / 29 |
| sr | self (n=40000) | random, strong hash | 18.4 | -0.1% [-0.22, 0.09] | -0.0% [-0.90, 0.81] | 100.00% | 100.00% | -0.00 [-0.008, 0.000] | 14 / 30 |
| sr | self (n=40000) | miniception | 15.7 | -14.5% [-14.66, -14.40] | -13.2% [-14.00, -12.47] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 15 / 31 |
| sr | self (n=40000) | double decycling | 16.4 | -11.2% [-11.30, -11.02] | -12.5% [-13.28, -11.77] | 100.00% | 100.00% | -0.01 [-0.013, 0.000] | 15 / 30 |
| sr | self (n=40000) | GreedyMini (min over strands) | 15.7 | -14.6% [-14.70, -14.46] | -14.3% [-15.05, -13.46] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 15 / 30 |
| sr | self (n=40000) | mod-minimizer | 15.3 | -16.9% [-17.05, -16.81] | -19.5% [-20.21, -18.72] | 100.00% | 99.98% | -0.02 [-0.033, -0.008] | 10 / 27 |
| sr | self (n=40000) | open-closed mod-minimizer | 14.5 | -21.3% [-21.44, -21.22] | -21.4% [-22.15, -20.72] | 100.00% | 99.98% | -0.02 [-0.030, -0.005] | 10 / 27 |
| sr | self (n=40000) | open-closed minimizer | 14.8 | -19.7% [-19.76, -19.53] | -23.4% [-24.06, -22.75] | 100.00% | 99.99% | -0.01 [-0.018, 0.000] | 15 / 31 |
| sr | self (n=40000) | open syncmer | 17.2 | -6.8% [-6.90, -6.64] | -8.1% [-8.90, -7.37] | 100.00% | 100.00% | -0.00 [-0.008, 0.000] | 8 / 25 |
| sr | self (n=40000) | closed syncmer | 20.2 | +9.4% [9.26, 9.59] | +11.3% [10.27, 12.31] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 9 / 25 |
| map-ont | self (n=800) | minimap2 2.31 (baseline) | 4326.0 | — | — | 100.00% | 100.00% | — | 162 / 171 |
| map-ont | self (n=800) | minimap2 order, width-adapted invertible hash | 4201.4 | -2.9% [-2.98, -2.78] | -6.4% [-7.55, -5.24] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 163 / 170 |
| map-ont | self (n=800) | random, strong hash | 4208.7 | -2.7% [-2.80, -2.62] | -2.5% [-3.18, -1.81] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 165 / 173 |
| map-ont | self (n=800) | miniception | 3657.2 | -15.5% [-15.55, -15.37] | -8.8% [-9.59, -7.98] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 172 / 181 |
| map-ont | self (n=800) | double decycling | 3674.0 | -15.1% [-15.15, -15.00] | -13.3% [-13.83, -12.70] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 168 / 172 |
| map-ont | self (n=800) | GreedyMini (min over strands) | 3885.2 | -10.2% [-10.30, -10.08] | -6.4% [-7.05, -5.72] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 167 / 177 |
| map-ont | self (n=800) | mod-minimizer | 3392.2 | -21.6% [-21.67, -21.50] | -19.3% [-20.08, -18.62] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 170 / 189 |
| map-ont | self (n=800) | open-closed mod-minimizer | 3303.4 | -23.6% [-23.74, -23.54] | -23.7% [-24.91, -22.53] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 167 / 184 |
| map-ont | self (n=800) | open-closed minimizer | 3440.3 | -20.5% [-20.56, -20.40] | -15.5% [-16.08, -14.98] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 164 / 182 |
| map-ont | self (n=800) | open syncmer | 4676.5 | +8.1% [7.99, 8.23] | +16.1% [15.37, 16.75] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 161 / 182 |
| map-ont | self (n=800) | closed syncmer | 4406.6 | +1.9% [1.73, 2.00] | +21.4% [20.46, 22.31] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 169 / 185 |
| map-ont | pbself (n=717) | minimap2 2.31 (baseline) | 4401.1 | — | — | 100.00% | 100.00% | — | 199 / 444 |
| map-ont | pbself (n=717) | minimap2 order, width-adapted invertible hash | 4279.0 | -2.8% [-2.88, -2.68] | -5.2% [-6.92, -4.10] | 100.00% | 99.86% | -0.14 [-0.418, 0.000] | 196 / 448 |
| map-ont | pbself (n=717) | random, strong hash | 4278.9 | -2.8% [-2.87, -2.67] | -3.2% [-3.86, -2.43] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 205 / 423 |
| map-ont | pbself (n=717) | miniception | 3722.9 | -15.4% [-15.53, -15.28] | -8.9% [-9.65, -8.41] | 100.00% | 99.86% | -0.14 [-0.418, 0.000] | 211 / 442 |
| map-ont | pbself (n=717) | double decycling | 3733.3 | -15.2% [-15.27, -15.08] | -13.1% [-13.73, -12.43] | 100.00% | 99.72% | -0.28 [-0.697, 0.000] | 208 / 436 |
| map-ont | pbself (n=717) | GreedyMini (min over strands) | 3946.5 | -10.3% [-10.46, -10.20] | -6.7% [-7.47, -6.00] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 201 / 433 |
| map-ont | pbself (n=717) | mod-minimizer | 3452.2 | -21.6% [-21.68, -21.45] | -18.9% [-19.79, -18.19] | 100.00% | 99.86% | -0.14 [-0.418, 0.000] | 207 / 445 |
| map-ont | pbself (n=717) | open-closed mod-minimizer | 3359.6 | -23.7% [-23.80, -23.54] | -22.8% [-24.32, -21.66] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 208 / 446 |
| map-ont | pbself (n=717) | open-closed minimizer | 3498.1 | -20.5% [-20.63, -20.41] | -15.8% [-16.26, -15.33] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 196 / 450 |
| map-ont | pbself (n=717) | open syncmer | 4758.5 | +8.1% [7.93, 8.32] | +16.1% [15.42, 16.93] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 189 / 381 |
| map-ont | pbself (n=717) | closed syncmer | 4485.5 | +1.9% [1.72, 2.12] | +21.4% [20.49, 22.27] | 100.00% | 99.86% | -0.14 [-0.418, 0.000] | 205 / 411 |
| map-hifi | self (n=1500) | minimap2 2.31 (baseline) | 1444.1 | — | — | 100.00% | 100.00% | — | 113 / 134 |
| map-hifi | self (n=1500) | random, strong hash | 1446.5 | +0.2% [0.05, 0.27] | +0.6% [-0.22, 1.52] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 115 / 135 |
| map-hifi | self (n=1500) | miniception | 1309.2 | -9.3% [-9.43, -9.25] | -6.6% [-7.37, -5.82] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 119 / 145 |
| map-hifi | self (n=1500) | double decycling | 1255.7 | -13.1% [-13.13, -12.97] | -8.9% [-9.64, -8.14] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 114 / 145 |
| map-hifi | self (n=1500) | GreedyMini (min over strands) | 1270.5 | -12.0% [-12.12, -11.92] | -5.1% [-6.52, -3.84] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 118 / 150 |
| map-hifi | self (n=1500) | mod-minimizer | 1447.0 | +0.2% [0.08, 0.32] | -0.0% [-0.88, 0.95] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 116 / 135 |
| map-hifi | self (n=1500) | open-closed mod-minimizer | 1170.6 | -18.9% [-19.03, -18.86] | -14.9% [-15.86, -14.12] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 123 / 155 |
| map-hifi | self (n=1500) | open-closed minimizer | 1170.1 | -19.0% [-19.06, -18.90] | -14.3% [-15.23, -13.48] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 124 / 155 |
| map-hifi | self (n=1500) | open syncmer | 1604.4 | +11.1% [10.98, 11.21] | +20.9% [19.87, 21.93] | 100.00% | 100.00% | +0.00 [0.000, 0.000] | 114 / 145 |

**Table 4. Unique vs repeat strata** (NA19240 reads and tiles; repeat = ≥ 50 % of the baseline's seeds above mid_occ).

| preset | stratum | scheme | n | Δ true (95% CI) | ≥1 true | Δ≥1 pp (95% CI) | Δ≥3 pp (95% CI) | longest gap p95 (bp) |
|---|---|---|---|---|---|---|---|---|
| sr | hap:uniq | minimap2 2.31 (baseline) | 39084 | — | 99.98% | — | — | 35 |
| sr | hap:uniq | random, strong hash | 39084 | +0.1% [-0.07, 0.25] | 99.97% | -0.01 [-0.021, 0.005] | -0.02 [-0.054, 0.015] | 35 |
| sr | hap:uniq | mod-minimizer | 39084 | -16.7% [-16.84, -16.57] | 99.98% | -0.00 [-0.013, 0.008] | -0.05 [-0.092, -0.018] | 32 |
| sr | hap:uniq | open-closed mod-minimizer | 39084 | -21.1% [-21.23, -20.99] | 99.97% | -0.01 [-0.020, 0.005] | -0.11 [-0.154, -0.067] | 32 |
| sr | hap:uniq | open-closed minimizer | 39084 | -19.3% [-19.46, -19.22] | 99.97% | -0.01 [-0.018, 0.003] | -0.09 [-0.125, -0.051] | 35 |
| sr | hap:uniq | closed syncmer | 39084 | +9.7% [9.54, 9.90] | 99.97% | -0.01 [-0.020, 0.005] | +0.01 [-0.018, 0.043] | 30 |
| sr | hap:uniq | open syncmer | 39084 | -6.3% [-6.40, -6.10] | 99.99% | +0.01 [-0.008, 0.018] | -0.03 [-0.064, 0.003] | 30 |
| sr | hap:rep | minimap2 2.31 (baseline) | 254 | — | 98.03% | — | — | 92 |
| sr | hap:rep | random, strong hash | 254 | -9.6% [-13.25, -5.98] | 98.43% | +0.39 [0.000, 1.181] | -2.36 [-6.693, 1.575] | 99 |
| sr | hap:rep | mod-minimizer | 254 | -26.6% [-30.17, -22.89] | 96.06% | -1.97 [-4.331, 0.010] | -14.96 [-20.079, -9.843] | 121 |
| sr | hap:rep | open-closed mod-minimizer | 254 | -28.7% [-31.59, -25.82] | 96.85% | -1.18 [-3.150, 0.787] | -12.20 [-17.323, -6.693] | 114 |
| sr | hap:rep | open-closed minimizer | 254 | -30.3% [-33.22, -27.33] | 98.82% | +0.79 [0.000, 1.969] | -7.09 [-10.630, -3.543] | 95 |
| sr | hap:rep | closed syncmer | 254 | +12.0% [7.89, 16.12] | 98.82% | +0.79 [0.000, 1.969] | +1.57 [-1.969, 5.118] | 94 |
| sr | hap:rep | open syncmer | 254 | -11.5% [-14.93, -8.05] | 99.21% | +1.18 [0.000, 2.756] | -0.79 [-4.331, 2.372] | 90 |
| map-ont | hap:uniq | minimap2 2.31 (baseline) | 755 | — | 100.00% | — | — | 309 |
| map-ont | hap:uniq | random, strong hash | 755 | -2.9% [-2.93, -2.78] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 318 |
| map-ont | hap:uniq | mod-minimizer | 755 | -21.5% [-21.56, -21.41] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 319 |
| map-ont | hap:uniq | open-closed mod-minimizer | 755 | -23.5% [-23.56, -23.39] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 318 |
| map-ont | hap:uniq | open-closed minimizer | 755 | -20.4% [-20.45, -20.30] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 315 |
| map-ont | hap:uniq | closed syncmer | 755 | +1.9% [1.82, 2.05] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 310 |
| map-ont | hap:uniq | open syncmer | 755 | +8.4% [8.25, 8.49] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 306 |
| map-ont | hap:rep | minimap2 2.31 (baseline) | 40 | — | 95.00% | — | — | 31147 |
| map-ont | hap:rep | random, strong hash | 40 | +0.2% [-5.77, 2.96] | 95.00% | +0.00 [0.000, 0.000] | -7.50 [-15.000, 0.000] | 23272 |
| map-ont | hap:rep | mod-minimizer | 40 | -17.0% [-21.60, -14.03] | 95.00% | +0.00 [0.000, 0.000] | -15.00 [-27.500, -5.000] | 31147 |
| map-ont | hap:rep | open-closed mod-minimizer | 40 | -34.8% [-40.98, -31.65] | 75.00% | -20.00 [-32.500, -7.500] | -42.50 [-57.500, -27.500] | 31319 |
| map-ont | hap:rep | open-closed minimizer | 40 | -16.3% [-21.61, -13.28] | 95.00% | +0.00 [0.000, 0.000] | -15.00 [-27.500, -5.000] | 31147 |
| map-ont | hap:rep | closed syncmer | 40 | +5.2% [-2.34, 13.87] | 95.00% | +0.00 [0.000, 0.000] | -12.50 [-25.000, -5.000] | 31147 |
| map-ont | hap:rep | open syncmer | 40 | +2.7% [-1.45, 6.39] | 95.00% | +0.00 [0.000, 0.000] | -2.50 [-7.500, 0.000] | 31147 |
| map-hifi | hap:uniq | minimap2 2.31 (baseline) | 1402 | — | 100.00% | — | — | 157 |
| map-hifi | hap:uniq | random, strong hash | 1402 | +0.1% [0.00, 0.22] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 149 |
| map-hifi | hap:uniq | mod-minimizer | 1402 | +0.2% [0.06, 0.26] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 149 |
| map-hifi | hap:uniq | open-closed mod-minimizer | 1402 | -18.6% [-18.64, -18.50] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 154 |
| map-hifi | hap:uniq | open-closed minimizer | 1402 | -18.6% [-18.69, -18.56] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 154 |
| map-hifi | hap:uniq | open syncmer | 1402 | +11.7% [11.64, 11.85] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 148 |
| map-hifi | hap:rep | minimap2 2.31 (baseline) | 80 | — | 90.00% | — | — | 15628 |
| map-hifi | hap:rep | random, strong hash | 80 | +0.1% [-4.83, 3.88] | 80.00% | -10.00 [-17.500, -2.500] | -8.75 [-16.250, -2.500] | 18438 |
| map-hifi | hap:rep | mod-minimizer | 80 | -2.8% [-8.01, 1.47] | 80.00% | -10.00 [-17.500, -2.500] | -8.75 [-15.000, -2.500] | 18438 |
| map-hifi | hap:rep | open-closed mod-minimizer | 80 | -19.7% [-25.84, -15.45] | 73.75% | -16.25 [-25.000, -8.750] | -30.00 [-41.250, -20.000] | 19575 |
| map-hifi | hap:rep | open-closed minimizer | 80 | -16.7% [-23.14, -12.75] | 73.75% | -16.25 [-25.000, -8.750] | -30.00 [-40.000, -20.000] | 19575 |
| map-hifi | hap:rep | open syncmer | 80 | +18.0% [7.79, 24.44] | 72.50% | -17.50 [-26.250, -10.000] | -26.25 [-36.250, -16.250] | 19575 |
| map-hifi | asm:uniq | minimap2 2.31 (baseline) | 3124 | — | 100.00% | — | — | 158 |
| map-hifi | asm:uniq | random, strong hash | 3124 | +0.1% [-0.00, 0.12] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 154 |
| map-hifi | asm:uniq | mod-minimizer | 3124 | +0.1% [0.03, 0.15] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 154 |
| map-hifi | asm:uniq | open-closed mod-minimizer | 3124 | -18.8% [-18.85, -18.76] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 174 |
| map-hifi | asm:uniq | open-closed minimizer | 3124 | -18.8% [-18.88, -18.80] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 174 |
| map-hifi | asm:uniq | open syncmer | 3124 | +11.3% [11.26, 11.38] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 153 |
| map-hifi | asm:rep | minimap2 2.31 (baseline) | 191 | — | 87.96% | — | — | 20000 |
| map-hifi | asm:rep | random, strong hash | 191 | +1.0% [-2.37, 4.10] | 81.68% | -6.28 [-10.995, -2.618] | -7.85 [-12.042, -3.665] | 20000 |
| map-hifi | asm:rep | mod-minimizer | 191 | -1.3% [-4.68, 1.80] | 81.15% | -6.81 [-11.008, -2.618] | -7.85 [-12.042, -3.665] | 20000 |
| map-hifi | asm:rep | open-closed mod-minimizer | 191 | -20.6% [-24.49, -17.33] | 75.39% | -12.57 [-17.801, -8.377] | -20.94 [-26.702, -15.183] | 20000 |
| map-hifi | asm:rep | open-closed minimizer | 191 | -18.3% [-22.32, -14.69] | 75.39% | -12.57 [-17.291, -7.853] | -20.94 [-26.702, -15.183] | 20000 |
| map-hifi | asm:rep | open syncmer | 191 | +15.0% [8.30, 20.38] | 75.39% | -12.57 [-17.277, -7.853] | -17.80 [-23.560, -12.042] | 20000 |
| asm20 | asm:uniq | minimap2 2.31 (baseline) | 3124 | — | 100.00% | — | — | 148 |
| asm20 | asm:uniq | random, strong hash | 3124 | +0.1% [0.05, 0.13] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 145 |
| asm20 | asm:uniq | mod-minimizer | 3124 | -21.3% [-21.31, -21.25] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 156 |
| asm20 | asm:uniq | open-closed mod-minimizer | 3124 | -24.7% [-24.76, -24.70] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 153 |
| asm20 | asm:uniq | open-closed minimizer | 3124 | -18.9% [-18.89, -18.83] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 151 |
| asm20 | asm:uniq | closed syncmer | 3124 | +0.5% [0.43, 0.52] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 148 |
| asm20 | asm:uniq | open syncmer | 3124 | +9.9% [9.81, 9.90] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 148 |
| asm20 | asm:rep | minimap2 2.31 (baseline) | 191 | — | 80.63% | — | — | 20000 |
| asm20 | asm:rep | random, strong hash | 191 | -2.0% [-4.00, -0.03] | 79.06% | -1.57 [-4.712, 1.571] | +0.00 [-2.618, 3.141] | 20000 |
| asm20 | asm:rep | mod-minimizer | 191 | -15.9% [-17.76, -13.84] | 76.96% | -3.66 [-6.283, -1.047] | +0.00 [-2.618, 2.618] | 20000 |
| asm20 | asm:rep | open-closed mod-minimizer | 191 | -25.2% [-27.10, -23.44] | 76.96% | -3.66 [-7.330, -0.524] | -3.66 [-7.330, 0.000] | 20000 |
| asm20 | asm:rep | open-closed minimizer | 191 | -19.8% [-22.25, -16.99] | 75.92% | -4.71 [-7.853, -2.094] | -4.19 [-7.330, -1.047] | 20000 |
| asm20 | asm:rep | closed syncmer | 191 | +9.7% [6.62, 14.03] | 76.96% | -3.66 [-6.283, -1.558] | +0.00 [-2.618, 3.141] | 20000 |
| asm20 | asm:rep | open syncmer | 191 | +0.5% [-2.06, 3.01] | 75.92% | -4.71 [-7.853, -2.094] | -1.57 [-4.712, 1.571] | 20000 |
| lr:hqae | asm:uniq | minimap2 2.31 (baseline) | 3126 | — | 100.00% | — | — | 226 |
| lr:hqae | asm:uniq | random, strong hash | 3126 | -0.3% [-0.37, -0.15] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 246 |
| lr:hqae | asm:uniq | mod-minimizer | 3126 | -0.2% [-0.26, -0.06] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 246 |
| lr:hqae | asm:uniq | open-closed mod-minimizer | 3126 | -12.9% [-12.96, -12.80] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 238 |
| lr:hqae | asm:uniq | open-closed minimizer | 3126 | -13.0% [-13.07, -12.91] | 100.00% | +0.00 [0.000, 0.000] | +0.00 [0.000, 0.000] | 234 |
| lr:hqae | asm:rep | minimap2 2.31 (baseline) | 189 | — | 75.13% | — | — | 20000 |
| lr:hqae | asm:rep | random, strong hash | 189 | +30.9% [20.80, 44.80] | 78.84% | +3.70 [-2.116, 8.995] | +14.29 [8.466, 20.106] | 20000 |
| lr:hqae | asm:rep | mod-minimizer | 189 | +31.1% [19.42, 44.98] | 79.37% | +4.23 [-1.058, 9.524] | +14.29 [8.466, 20.106] | 20000 |
| lr:hqae | asm:rep | open-closed mod-minimizer | 189 | -9.0% [-17.81, -1.04] | 70.37% | -4.76 [-8.995, -0.529] | -4.76 [-9.524, -0.529] | 20000 |
| lr:hqae | asm:rep | open-closed minimizer | 189 | -9.1% [-17.93, -0.85] | 70.37% | -4.76 [-9.008, -0.529] | -4.76 [-8.995, 0.000] | 20000 |

**Table 5. Stock minimap2 2.31 vs 2.31 + `--mod-oc`** (CHM13 chr20 index, 1 thread, paftools.js mapeval; err = fraction of mapped reads not overlapping the true locus).

| preset | variant | reads | .mmi bytes | index CPU s | index RSS GB | map CPU s | map RSS GB | MAPQ60 mapped (err) | MAPQ≥1 mapped (err) | all mapped (err) |
|---|---|---|---|---|---|---|---|---|---|---|
| sr | stock 2.31 | 79281 | 197,061,282 | 2.33 | 0.911 | 5.23 | 0.349 | 72363 (0.001%) | 74860 (0.319%) | 79281 (5.371%) |
| sr | 2.31 + --mod-oc | 79281 | 153,923,482 | 3.37 | 0.745 | 4.18 | 0.313 | 72585 (0.001%) | 74832 (0.271%) | 79281 (5.377%) |
| map-ont | stock 2.31 | 1574 | 210,102,234 | 2.50 | 0.969 | 37.66 | 0.475 | 1517 (0.000%) | 1546 (0.259%) | 1574 (1.779%) |
| map-ont | 2.31 + --mod-oc | 1574 | 166,030,442 | 3.25 | 0.782 | 26.33 | 0.327 | 1515 (0.000%) | 1547 (0.259%) | 1574 (1.779%) |
| map-hifi | stock 2.31 | 2942 | 130,955,730 | 1.66 | 0.533 | 14.92 | 0.254 | 2826 (0.000%) | 2902 (0.345%) | 2942 (1.428%) |
| map-hifi | 2.31 + --mod-oc | 2942 | 112,287,554 | 3.07 | 0.499 | 13.07 | 0.248 | 2826 (0.000%) | 2895 (0.138%) | 2942 (1.462%) |

## Limitations and open issues
* The seed metrics use minimap2's hard mid_occ cut-off. `mm_seed_select`'s occ_dist rescue of high-occurrence seeds is not modelled. The repeat-stratum losses (0.6–5.8 % of NA19240 reads) are therefore an upper bound on harm. The aligner-level test, which includes the rescue, showed no accuracy loss on chr20 or genome-wide for HiFi.
* All reads are simulated. The ONT Q20 model's error split is illustrative. pbsim3's ERRHMM-ONT-HQ at `--accuracy-mean 0.99` realised 2.95 % errors. HiFi follows the Wenger et al. 2019 composition, not a modern Revio/DeepConsensus profile.
* Truth for NA19240 reads comes from a minimap2 asm5 alignment of each haplotype to CHM13. About 1 % of bases are unliftable, and reads with < 50 % liftable seeds (1.7 % of Illumina hap reads) are excluded.
* GreedyMini orders exist only for k ≤ 16 and binary orders. k = 19/21 used the documented k′-prefix extension and lr:hqae had none. Canonicalization used min over strands, not GreedyCan-generated canonical orders, so GreedyMini is under-represented here.
* The `--mod-oc` sketch routine is scalar and unoptimised: index construction is 1.3–1.9× slower. A SIMD implementation in the style of SimdMinimizers ([Groot Koerkamp & Martayan 2025, preprint](https://doi.org/10.1101/2025.01.27.634998)) is the obvious next step.
* The new .mmi magic is a format decision for the maintainer. Draft PR lh3/minimap2#1361 documents the silent-misuse hazard that motivates it.
* map-pb (HPC) and splice presets were not studied, beyond splice's hash-quality effect.

## Proposals
1. **Open-closed mod-minimizer as a persisted minimap2 index option** (index/sketching; prototype; `minimap2_sampling.patch`). The closest prior work is draft PR lh3/minimap2#1361, a sketch-time plain mod-minimizer with no gain when t = k and no persisted index. Litmus: L1-A above.
2. **Width-adapted invertible hash for k ≤ 15** (index/sketching; open). Litmus: L1-B.
3. **Canonical-by-construction low-density orders** (index/sketching; open). Make GreedyCan canonical orders and a parity-aware decycling variant for (15, 10), (19, 19) and (21, 11). Litmus: canonical chr20 density at least 5 % below modoc (0.1418 / 0.0816 / 0.1235) with strand symmetry ≥ 0.999, then L1-A.
4. **Equal-density syncmer or mod sampling for short reads** (short-read alignment; implemented in strobealign, not adopted in minimap2). Litmus: L1-C.
5. **Density-aware high-occurrence handling** (cross-cutting; open). Scale occ_dist and mid_occ with density so repeat reads keep anchors. Litmus: for NA19240 repeat-stratum reads (Table 4), the fraction with ≥ 1 retained true seed under modoc plus occ_dist rescue is within 2 pp of stock minimap2 with rescue (stock without rescue: 90.0 % HiFi, 95.0 % ONT).

## Files and reproduction

Run `bash L1_pipeline.sh WORK_DIR [--wholegenome]`. Single-threaded; Python steps peaked at 3–6 GB RSS, and the whole-genome step at 8.7 GB. Inputs are downloaded by `L1_fetch_data.sh`, with md5s written to `data/MD5SUMS.txt`: CHM13v2.0 md5 06caded0055ec647a69acb708e13beff, GRCh38 chr20 .fa.gz md5 1d42fbb5014b165373aca053f4f24061, minimap2-2.31.tar.gz md5 ca5ddc443d95f040637fb7bc3ff17bb1. Tool versions: minimap2 2.31-r1302, samtools/htslib 1.25, pbsim3 3.0.5, k8 1.2, Python 3.11.18, numba 0.68.1, numpy 2.4.6, pandas 3.1.1, matplotlib 3.11.2, Apple clang 21.0.0 (macOS arm64).

| file | content |
|---|---|
| [L1_schemes.py](L1_schemes.py) | sampling schemes (numba): exact mm_sketch port, canonical miniception / double decycling / GreedyMini / mod / modoc / OC / syncmers, hash_w |
| [L1_simulate.py](L1_simulate.py) | read simulator with per-base truth (Illumina, HiFi, ONT models; pbsim3 import) |
| [L1_eval.py](L1_eval.py) | index construction (minimap2 mid_occ rule), liftover, per-read true-seed evaluation |
| [L1_run.py](L1_run.py) | experiment driver (parts main, grch38, blocks, curves, hashw) |
| [L1_stats.py](L1_stats.py) | paired / block bootstrap statistics → L1_results.tsv |
| [L1_tables.py](L1_tables.py) | markdown tables |
| [L1_figure.py](L1_figure.py) | figure |
| [L1_validate.py](L1_validate.py) | validation suite → L1_validation.tsv |
| [L1_make_fastq.py](L1_make_fastq.py) | mapeval FASTQs (truth in read names) |
| [L1_patch_eval.sh](L1_patch_eval.sh) | stock vs patched minimap2 on chr20 (index, timing, mapeval) |
| [L1_wholegenome.sh](L1_wholegenome.sh) | whole-genome CHM13v2.0 index/mapping comparison |
| [L1_fetch_data.sh](L1_fetch_data.sh) | data download + liftover alignments |
| [L1_pipeline.sh](L1_pipeline.sh) | end-to-end driver |
| [L1_sketch_harness.c](L1_sketch_harness.c) | C harness around minimap2 sketch.c |
| [minimap2_sampling.patch](minimap2_sampling.patch) | minimap2 2.31 patch: --mod-oc (flag 0x8, magic MMI\3), 241 lines |
| [L1_results.tsv](L1_results.tsv) | all index and per-read metrics with CIs (357 rows; groups incl. :uniq/:rep strata) |
| [L1_index_chm13.tsv](L1_index_chm13.tsv) | index statistics, CHM13 chr20 |
| [L1_index_grch38.tsv](L1_index_grch38.tsv) | index statistics, GRCh38 chr20 |
| [L1_curves.tsv](L1_curves.tsv) | equal-density curves |
| [L1_equal_density.tsv](L1_equal_density.tsv) | curves interpolated at the preset density |
| [L1_validation.tsv](L1_validation.tsv) | validation results |
| [L1_patch_eval.tsv](L1_patch_eval.tsv) | patched vs stock, chr20 |
| [L1_wholegenome.tsv](L1_wholegenome.tsv) | patched vs stock, whole genome |
| [L1_mapeval.txt](L1_mapeval.txt) | raw paftools.js mapeval output |
| [L1_readsets.json](L1_readsets.json) | read-set sizes and calibrated HP constants |
| [L1_references.json](../../references/tracks/L1_references.json) | verified bibliography (DOIs) |
| [L1_tables.md](L1_tables.md) | rendered tables |

## References

- Roberts; Hayes; Hunt et al. (2004). Reducing storage requirements for biological sequence comparison. Bioinformatics. https://doi.org/10.1093/bioinformatics/bth408
- Schleimer; Wilkerson; Aiken (2003). Winnowing. Proceedings of the 2003 ACM SIGMOD international conference on Management of data. https://doi.org/10.1145/872757.872770
- Li (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty191
- Li (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics. https://doi.org/10.1093/bioinformatics/btab705
- Marçais; Pellow; Bork et al. (2017). Improving the performance of minimizers and winnowing schemes. Bioinformatics. https://doi.org/10.1093/bioinformatics/btx235
- Marçais; DeBlasio; Kingsford (2018). Asymptotically optimal minimizers schemes. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty258
- Zheng; Kingsford; Marçais (2020). Improved design and analysis of practical minimizers. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa472
- Pellow; Pu; Ekim et al. (2023). Efficient minimizer orders for large values of k using minimum decycling sets. Genome Research. https://doi.org/10.1101/gr.277644.123
- Groot Koerkamp; Pibiri (2024). The mod-minimizer: a simple and efficient sampling algorithm for long k -mers. bioRxiv (preprint). https://doi.org/10.1101/2024.05.25.595898
- Groot Koerkamp; Liu; Pibiri (2025). The open-closed mod-minimizer algorithm. Algorithms for Molecular Biology. https://doi.org/10.1186/s13015-025-00270-0
- Kille; Koerkamp; McAdams et al. (2024). A near-tight lower bound on the density of forward sampling schemes. bioRxiv (preprint). https://doi.org/10.1101/2024.09.06.611668
- Golan; Tziony; Kraus et al. (2025). GreedyMini: generating low-density DNA minimizers. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaf251
- Shur; Tziony; Orenstein (2026). Generating Minimum-Density Minimizers. Lecture Notes in Computer Science. https://doi.org/10.1007/978-3-032-42866-0_20
- Shur; Tziony; Orenstein (2026). 10-minimizers: a promising class of constant-space minimizers. bioRxiv (preprint). https://doi.org/10.64898/2026.03.16.712052
- Edgar (2021). Syncmers are more sensitive than minimizers for selecting conserved k ‑mers in biological sequences. PeerJ. https://doi.org/10.7717/peerj.10805
- Shaw; Yu (2021). Theory of local k-mer selection with applications to long-read alignment. Bioinformatics. https://doi.org/10.1093/bioinformatics/btab790
- Dutta; Pellow; Shamir (2022). Parameterized syncmer schemes improve long-read mapping. PLOS Computational Biology. https://doi.org/10.1371/journal.pcbi.1010638
- Sahlin (2021). Effective sequence similarity detection with strobemers. Genome Research. https://doi.org/10.1101/gr.275648.121
- Sahlin (2022). Strobealign: flexible seed size enables ultra-fast and accurate read alignment. Genome Biology. https://doi.org/10.1186/s13059-022-02831-7
- Jain; Rhie; Zhang et al. (2020). Weighted minimizer sampling improves long read mapping. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa435
- Jain; Rhie; Hansen et al. (2022). Long-read mapping to repetitive reference sequences using Winnowmap2. Nature Methods. https://doi.org/10.1038/s41592-022-01457-8
- Shaw; Yu (2023). Proving sequence aligners can guarantee accuracy in almost O ( m log n ) time through an average-case analysis of the seed-chain-extend heuristic. Genome Research. https://doi.org/10.1101/gr.277637.122
- Zheng; Marçais; Kingsford (2023). Creating and Using Minimizer Sketches in Computational Genomics. Journal of Computational Biology. https://doi.org/10.1089/cmb.2023.0094
- Groot Koerkamp; Martayan (2025). SimdMinimizers: Computing random minimizers, fast. bioRxiv (preprint). https://doi.org/10.1101/2025.01.27.634998
- Nurk; Koren; Rhie et al. (2022). The complete sequence of a human genome. Science. https://doi.org/10.1126/science.abj6987
- Liao; Asri; Ebler et al. (2023). A draft human pangenome reference. Nature. https://doi.org/10.1038/s41586-023-05896-x
- Logsdon; Ebert; Audano et al. (2025). Complex genetic variation in nearly complete human genomes. Nature. https://doi.org/10.1038/s41586-025-09140-6
- Gustafson; Gibson; Damaraju et al. (2024). High-coverage nanopore sequencing of samples from the 1000 Genomes Project to build a comprehensive catalog of human genetic variation. Genome Research. https://doi.org/10.1101/gr.279273.124
- Wenger; Peluso; Rowell et al. (2019). Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. Nature Biotechnology. https://doi.org/10.1038/s41587-019-0217-9
- Sereika; Kirkegaard; Karst et al. (2022). Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. Nature Methods. https://doi.org/10.1038/s41592-022-01539-7
- Ono; Hamada; Asai (2022). PBSIM3: a simulator for all types of PacBio and ONT long reads. NAR Genomics and Bioinformatics. https://doi.org/10.1093/nargab/lqac092
- WABI 2024 proceedings version of the mod-minimizer paper (LIPIcs). https://doi.org/10.4230/LIPIcs.WABI.2024.11

*Corrected on 2026-11-12 after an independent fact-check (`factcheck/factcheck_report.md`).*
