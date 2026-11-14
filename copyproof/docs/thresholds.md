# Thresholds: what each value is and where it came from

Every numeric choice in copyproof, sorted by how it was arrived at. "Chosen" means a value was
picked on judgement and kept because results looked reasonable; it was not optimised and has no
derivation. Where the notebook or an experiment gives the reason, the entry is named.

None of these was fixed before the data were seen, and several were changed on the same data the
results are reported on (rungs 1 to 4). A second sample with every value frozen beforehand is the
test still missing (GUIDE §8 rung 6).

The classifier's thresholds were varied one at a time in
[`experiments/2026-11-13_threshold_sensitivity`](../experiments/2026-11-13_threshold_sensitivity/README.md);
"plateau" below means no count there moved by more than a few.

## Derived, fitted or measured

| What | Value | Where it comes from |
|---|---|---|
| Depth per copy | per run (13.4 at 36x) | half the median read count of k-mers found once in each haplotype (`windows.depth_per_copy`) |
| GC curve | per run | median depth of two-copy k-mers per 1% GC bin of 1-kb windows |
| Extra-Poisson noise per class (`omega`) | per run (1.0 to 4.4) | fitted, after the covariance of neighbouring k-mers is taken out exactly (NOTEBOOK d) |
| Regional wander of depth (`psi`) | per run (7.7% at 36x) | fitted from the two-copy class |
| k-mer starts per read (`reff`) | per run (120 for 150-bp reads) | the reads' mean length in the scan, minus k, plus 1 |

## Set from a count on real data

In this group the rule exists because of a count; where the exact value was still a choice, the
row says so.

| Parameter | Value | The count behind it | Record |
|---|---|---|---|
| `absent_min_windows` | 20 | shorter stretches at under a quarter of their depth: 263 on the truth genome, 1 real | NOTEBOOK e |
| `use_baseline` | off | a local depth baseline took the truth genome's whole-copy calls from 17 to 40 | NOTEBOOK e |
| `fit_share` | 0.9 | with every class able to veto a fit, four large planted events were lost; sensitivity: at 1.0, 7 planted signals and 12 surplus copies fewer | NOTEBOOK j |
| `max_wander` | 2 sd | a +1 call on the truth genome needed a regional factor of 0.79; plateau | NOTEBOOK j |
| `deep_ratio`, `deep_windows` | 1.5-fold, 8 windows | the truth genome has no ambiguous segment that long and deep, and 105 short dips to half depth (`ambiguous_by_length.tsv`) | rung 3 |
| `high_windows` | 5 | whole-copy losses under 5 windows: 13 on the truth genome (long reads hold the sequence at 10), 18 in an ordinary assembly of which the truth confirms 1 | rungs 3 and 4 |
| Weak one-sided missing sequence is not an issue | under 6 reads or 5 tracked k-mers | true 9% to 10% of the time where truth is known, against 46% to 51% for strong ones and 96% for strong missing alleles | rung 4, item 9 |
| ... unless in neither parent | both parents at 0 | 46 such events in HG01943, the known de novo junction first among them; 5 on the truth genome | rung 1, trio |
| Strong one-sided missing sequence is never high-confidence | | a first rule raised it when a parent's reads held the k-mers; on the truth genome that holds for 208 of 211 strong clusters | trio |
| Hapmers: both parents readable (`parent_depth`) | 0.7 to 1.5 of the parent's usual depth on two-copy k-mers | the rule came from the count (111 of 143 false switches on the truth genome lie in such windows); the two bounds were chosen and not varied. Effect: 143 segments to 32, with 1.7% of windows set aside | trio |
| Hapmers: density of a switch (`min_density`) | a tenth of the sequence's own site density | with the row above, 32 to 24 segments and 5.6 to 0.4 Mb. Any value from 0.1 to 0.2 gives the same 24; at 0.2 one planted 200-kb swap, whose ratio is 0.20, is lost | trio |
| Unlinked join is firm only next to another issue | within 20 kb | 34 of 42 scaffold joins of one release have no linking long read | rung 1 |
| Long-read verdicts | carried: 80% of the unread k-mers in at least 2 reads; not carried: at most 20%, with at least 5 reads | not tuned; checked afterwards: 115 of 119 not-carried runs are known issues of the truth genome, 74 of 957 carried ones | rung 3 |

## Chosen, with a sensitivity check

| Parameter | Value | What the check showed |
|---|---|---|
| `psi_local`: how far one k-mer class may sit from another in a segment | 4% | 3% is the same; at 6% a quarter of the collapses found in the assembly are lost. For scale: copies per class come out within 2.5% of the multiplicity genome-wide (NOTEBOOK d) |
| `max_z`: a class fits within this many sd | 3.5 | plateau between 3 and 4 |
| `margin`: log-likelihood to prefer one explanation | 3 nats (about 20 to 1) | **2 looks better on these data**: 13 more confirmed high-confidence calls on the assembly, 11 more on the planted genome, the truth genome unchanged. Left at 3 |
| `dominant_tol`: the commonest class lands this near a whole copy number | 0.35 | plateau between 0.25 and 0.45 |
| `min_share`: share of k-mers that must follow a whole-copy change | 0.5 | 0.4 the same; at 0.6, 5 collapses fewer |
| `min_seg_class`: k-mers a class needs to count in a segment | 20 | **10 looks better on these data**: ambiguous segments on the truth genome 143 to 49, more high-confidence calls on the other two. 40 is worse on all three. Left at 20 |
| `min_llr`: log-likelihood a segment needs | 20 nats | at 30: truth-genome issues 49 to 17, 11 planted signals lost, the assembly keeps 110 of 114 confirmed calls |

## Chosen, not varied

These need the pipeline to run again to vary, and were not.

| Where | Parameter | Value | Note |
|---|---|---|---|
| Engine | k | 31 | the pilots' junction k-mers and the usual hapmer size; other sizes were not tried (decision D2) |
| Engine | k-mers tracked | 1 in 4 | neighbouring k-mers share reads, so little is lost: a heterozygous site still makes about 8 tracked k-mers. `--sample 8` halves memory |
| Engine | read k-mers absent from the release are listed from | 3 reads | |
| Engine | high-copy k-mers (not used as evidence of missing sequence) | over 8 copies in the release, or read at over 20 copies' depth | |
| Windows | size | 1 kb | sets the smallest event: changes under 5 kb are mostly missed |
| Windows | counts winsorised per k-mer at | m + 4 copies | |
| Windows | highest multiplicity used (`cap`) | 8 | above it a window's k-mers are not used; over half of them above it masks the window |
| Model | noise | Student-t, 5 degrees of freedom | |
| Mask | usable k-mers per window; per class | 30; 4 | |
| Mask | GC factor; A/C/G/T share; shortest assessable island | 0.5; 0.5; 5 windows | the mask as a whole came from the first truth-genome runs (NOTEBOOK e); its values were not varied |
| HMM | prior that an event starts at a window; ends | 10^-6; 0.02 | |
| HMM | states | +1, -1, +2 copies for k-mers of multiplicity at least 1, 2, 3 or 5; factors 0, 0.25, 0.5, 0.75, 1.25, 1.5, 2 | changes of more than 2 copies are found as their nearest state and classified from the summed evidence |
| HMM | segments joined across | 2 windows | |
| Other copies | a copy holds at least; diagonal tolerance; unique k-mers needed | 20% of the shared k-mers and 20 of them; 2 kb; 20 | on the truth genome this step changed one segment |
| Classifier | "absent" | every class at a quarter of its depth or less | |
| Support | another channel's issue within; a join within | 1 kb; 5 kb | |
| S2 | solid; ceiling | half the depth per copy; 3.5 copies | |
| S2 | recruit a read with; skip it above | 2 missing k-mers; half of its k-mers high-copy | |
| S2 | anchors: occurrences in the release; co-linear k-mers per flank | at most 4; 3 | a flank that is repeated more than 4 times cannot anchor, which is why a junction at a high-copy repeat is found weakly (rung 1) |
| S2 | reads and tracked k-mers a cluster needs; each strand | 3 and 3; 1 read and 15% | the strand rule dropped 13,102 loci on the truth genome against 1,408 kept (NOTEBOOK f) |
| S2 | strong cluster | 6 reads and 5 tracked k-mers | |
| S2 | missing haplotype | 3 missing alleles within 10 kb; high confidence from 5 | |
| S3 | hapmer | in the child, at least 3 reads in one parent, none in the other | |
| S3 | wrong-parent sites assumed; switch prior | 2%; 10^-3 per window | |
| S3 | a switch needs; a sequence or piece is assessed from | 4 sites at 3 to 1; 20 sites | |
| S4 | join: bases fetched each side; a read must reach; reads needed near each side | 10 kb; 5 kb; 3 | |
| S4 | run: alignment must cover; clip counted from; within | 500 bp beyond both ends; 200 bp; 100 bp | |
