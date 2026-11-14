# How much do the results hang on the thresholds, and on the classifier's revision?

**Question.** Most of the classifier's numbers were chosen, not derived (`docs/thresholds.md`
says which). Do the results move when they do? And the classifier was revised during rung 4
(NOTEBOOK entry j): what did the revision change, measured the same way on all three kinds of
data?

**Inputs.** The cached segmentations of three finished runs: the truth genome (rung 3), the
planted genome (rung 2) and one ordinary assembly, `hprc_year2` (rung 4).

**Run.**

    scripts/dev/cache_segments.py RUN_DIR            # once per run: the HMM's segments
    scripts/dev/sensitivity.py --truth ... --planted ... --expected plant1.expected.bed \
        --asm .../hprc_year2 --asm-truth hprc_year2.eval.truth_cn_segments_raw.tsv --out sensitivity.tsv

The HMM is not run again; each setting re-classifies the same segments and re-applies the rules
that turn segments into issues. One threshold is changed at a time. "Before the revision" uses
the earlier classifier, kept in `scripts/dev/classifier_before_revision.py`, with today's issue
rules. Two pipeline steps are left out (the lookup of other copies, support from other channels),
so the numbers sit a little below those of the experiment records.

## Result

On the truth genome fewer is better. On the planted genome 138 signals are expected (duplications,
extra copies and deletions with their homologous stretches). On the assembly, 64 collapses and
169 surplus copies are assessable discordances by the truth genome's k-mers.

| Setting | Truth genome: whole-copy segments | ... ambiguous | ... issues, high / low | Planted: found (high) of 138 | Assembly: high-confidence confirmed | ... collapses found | ... surplus copies found |
|---|---|---|---|---|---|---|---|
| **before the revision** | 14 | 180 | 2 / 35 | 107 (81) | 89 of 91 | 8 of 64 | 123 of 169 |
| **as shipped** | 20 | 143 | 4 / 45 | 108 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `psi_local = 0.03` | 20 | 142 | 4 / 45 | 107 (82) | 116 of 123 | 19 of 64 | 131 of 169 |
| `psi_local = 0.06` | 18 | 151 | 3 / 43 | 106 (76) | 99 of 101 | 13 of 64 | 135 of 169 |
| `max_z = 3` | 18 | 134 | 4 / 40 | 106 (82) | 114 of 120 | 17 of 64 | 133 of 169 |
| `max_z = 4` | 17 | 150 | 4 / 46 | 107 (78) | 113 of 118 | 17 of 64 | 136 of 169 |
| `margin = 2` | 22 | 107 | 4 / 45 | 108 (92) | 127 of 134 | 19 of 64 | 135 of 169 |
| `margin = 4` | 18 | 163 | 3 / 37 | 104 (73) | 105 of 108 | 15 of 64 | 133 of 169 |
| `dominant_tol = 0.25` | 19 | 128 | 4 / 40 | 107 (81) | 112 of 117 | 16 of 64 | 132 of 169 |
| `dominant_tol = 0.45` | 20 | 145 | 4 / 45 | 108 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `fit_share = 0.8` | 20 | 180 | 4 / 60 | 108 (81) | 113 of 118 | 17 of 64 | 135 of 169 |
| `fit_share = 1` | 15 | 125 | 4 / 32 | 101 (77) | 106 of 111 | 17 of 64 | 123 of 169 |
| `max_wander = 1.5` | 20 | 147 | 4 / 45 | 108 (82) | 113 of 118 | 16 of 64 | 134 of 169 |
| `max_wander = 3` | 21 | 138 | 5 / 45 | 107 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `min_share = 0.4` | 20 | 143 | 4 / 45 | 108 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `min_share = 0.6` | 20 | 141 | 4 / 45 | 106 (80) | 110 of 115 | 12 of 64 | 134 of 169 |
| `min_seg_class = 10` | 17 | 49 | 4 / 22 | 108 (90) | 119 of 124 | 18 of 64 | 134 of 169 |
| `min_seg_class = 40` | 30 | 391 | 2 / 61 | 107 (77) | 102 of 106 | 15 of 64 | 134 of 169 |
| `min_llr = 30` | 9 | 29 | 2 / 15 | 97 (77) | 110 of 114 | 17 of 64 | 133 of 169 |
| `min_llr = 40` | 4 | 3 | 2 / 2 | 87 (70) | 108 of 112 | 15 of 64 | 133 of 169 |
| `deep_windows = 5` | 20 | 143 | 4 / 47 | 108 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `deep_windows = 12` | 20 | 143 | 4 / 45 | 96 (81) | 114 of 119 | 17 of 64 | 134 of 169 |
| `deep_ratio = 1.33` | 20 | 143 | 4 / 47 | 108 (81) | 114 of 119 | 17 of 64 | 135 of 169 |
| `deep_ratio = 2` | 20 | 143 | 4 / 45 | 91 (81) | 114 of 119 | 17 of 64 | 123 of 169 |
| `high_windows = 3` | 20 | 143 | 5 / 44 | 108 (82) | 115 of 121 | 17 of 64 | 135 of 169 |
| `high_windows = 8` | 20 | 143 | 1 / 48 | 108 (74) | 111 of 116 | 17 of 64 | 135 of 169 |

## What it says

1. **The revision bought recall on the real assembly at a small cost on the truth genome.**
   High-confidence issues confirmed on `hprc_year2` rose from 89 of 91 to 114 of 119; collapses
   found from 8 to 17 of 64; surplus copies from 123 to 135 of 169. The planted genome is
   unchanged (107 and 108 of 138). The truth genome went from 2 to 4 high-confidence issues and
   from 35 to 45 low-confidence ones.
2. **Most thresholds sit on a plateau.** Moving `max_z`, `dominant_tol`, `max_wander` or
   `min_share` (down) by a third changes no count by more than a few. `psi_local` at 3% is the
   same as at 4%; at 6% a quarter of the collapses are lost.
3. **Two rules that came from the data hold up.** With `fit_share` at 1 (every class may veto a
   fit, as before the revision) 7 planted signals and 12 surplus copies are lost; at 0.8 the
   truth genome gains 15 low-confidence issues. With `deep_ratio` at 2 or `deep_windows` at 12
   the planted genome loses 12 to 17 signals and nothing is gained on the truth genome.
4. **Two chosen values are not the best on these data, and were left alone.**
   - `min_seg_class` (k-mers a class needs to count in a segment's classification; 20). At 10,
     ambiguous segments on the truth genome fall from 143 to 49 and its low-confidence issues
     from 45 to 22, planted high-confidence calls rise from 81 to 90, and the assembly gains 5
     confirmed calls. At 40 everything is worse.
   - `margin` (3 nats). At 2, confirmed high-confidence calls on the assembly rise from 114 to
     127 and planted high-confidence calls from 81 to 92, with the truth genome unchanged at 4.

   Both would be changes made on the very data the results are reported on, after several
   rounds of that already. They are recorded here as the first things to try on a second sample
   with thresholds frozen beforehand (rung 6).
5. **`min_llr` trades small events for a quieter truth genome.** At 30 nats instead of 20 the
   truth genome's issues fall from 49 to 17 and the assembly keeps 110 of its 114 confirmed
   calls, but 11 planted signals are lost.

## Not covered

Thresholds of the HMM itself (the prior, the states, the window size), of the mask, and of the
missing-sequence, hapmer and long-read channels were not varied: they need the pipeline to run
again. `docs/thresholds.md` lists them with where each value came from.

## Files

- `sensitivity.tsv`: every count per setting.
