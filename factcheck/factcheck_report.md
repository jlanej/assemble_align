# Fact-check of the assemble_align review draft (first full draft), 12 November 2026

The draft's measured numbers are faithful to the pilot outputs. Every cell of Tables L1.1, L2.1, L5.1, L5.2 and L5.3 and every score in Table 3.1 matches its source table, and both corrections the report makes to the case repositories hold up under independent re-tests. The problems are scope and wording. The check produced 52 findings (10 errors, 21 mis-scoped, 9 unsupported, 2 stale, 10 style), and most can be fixed in one sentence. They are listed with exact quotes, evidence and replacement text in [factcheck_findings.tsv](factcheck_findings.tsv).

## Fix first

1. **L5 headline numbers lose their qualifiers (FC01-FC04).** The executive summary credits the calibrated recalibrator with keeping 28-58% of correct in-SD reads at ≤0.1% error. That figure comes from the raw GBM ranking, with the threshold set on held-out labels; the calibrated model keeps 0% for long reads. The 71-96% "paralog-better" share covers all wrong MAPQ-60 reads, not only those in SDs (80-97% in SDs). Table 2.1's "1.0-2.8% depending on the truth definition" is one truth definition; the strict ≥90%-colinear truth gives 0.12-0.33%.
2. **L1 headline claims (FC05-FC10).** "No loss of seeds in unique sequence" is wrong: true seeds per read fall 13-25%. What is unchanged is the share of long reads with ≥3 true seeds; Illumina loses 0.11 points. The CPU and wrong-mapping gains are chr20 results at MAPQ ≥ 1. Old binaries do not "refuse" the MMI\3 index; stock 2.31 parsed it as FASTA and mapped 0/400 reads with exit status 0. F6 also misdescribes PR #1361 (see the status table below). The Table L1.2 reads are about half CHM13-derived, not all NA19240.
3. **Minimizer lower bound (FC40-FC42).** Figure L1(a) plots ⌈(w+k)/w⌉/(w+k), and Part II quotes about 0.053 at w = k = 19 with "up to ~45%" avoidable. That simplified formula is not Kille et al.'s improved bound. Theorem 2 of that paper raises the bound to about 3/39 ≈ 0.077 at w = k = 19 (density factors 1.41-1.54 across the five presets). At most about 23% of entries are therefore avoidable at map-hifi, which fits L1's measured −18.8%.
4. **1kG-ONT naming and novelty (FC34, FC35).** Part I uses "1kG-ONT" for two cohorts. One is the Vienna set of 1,019 genomes at 16.9x, which has no whole-genome assemblies. The other is the 1kGP-LRSC set of about 500 genomes, which supplies the "12 HPRC-overlap samples". The LRSC bucket already holds hifiasm ONT-mode dual assemblies for 247 R10 samples (S3 listing; files dated 2025-09-15).
5. **Citation support (FC32, FC33).** Qtip is not "the only systematic learned predictor": L5's own background cites LoQuM and an AdaBoost-isotonic recalibrator. Locityper does not use a Li-Stephens copying model.
6. **Other errors.** F1 maps NYGC Illumina reads with the long-read preset `lr:hqae` (FC17). "24% more as gapped T2T scaffolds" should read 23%: 2,490 of 10,626 is 23.4% (FC21). Part II says NA19240's parents are in HPRC release 2; they are HGSVC3 samples (FC43).

## Case-repository corrections

| Correction | Repository statement | Verdict | Evidence |
|---|---|---|---|
| (a) HG01943 junction has a 2-bp `CT` insertion, not 1-bp microhomology | AssemblyReview `cases/hg01943/README.md`: "The junction has 1 bp of microhomology" | **Confirmed** (independent test) | The HPRC/Yale HiFi BAM was read over HTTPS. All 31 reads soft-clipped at chr2:48,057,041 continue `CT` + `GCTCTGGCAG…` (chr2:46,970,725 onward), and 33 reads clipped at 46,970,725 end `…TTCCACT`. Of 182 reads at the breakpoints, 32 carry the CT-junction 14-mer and none carry the plain join. The GRCh38 flanks (`tg` after the end, `ta` before the start) template neither base. The repository's "1 bp" matches the reference-level join, where chr2:48,057,041 and chr2:46,970,724 are both A. |
| (b) HG01103 trio paternal FASTA holds 16/16 junctions, not 5 | svreview- `examples/hg01103/README.md`: "It holds 5 of the 16 junctions (J6, J7, J8, J9, J10)" | **Confirmed**, with a wording fix | The hap1 FASTA was streamed: 217 contigs, 3.03 Gb. All 22-26 read-derived 31-mers of each of the 16 junctions sit on the contigs the track lists. In the published hap1 BAM, 142 of 217 contigs have records; the 75 absent contigs total 22.89 Mb, and no record has MAPQ < 5. Only 4 of the 11 hidden junctions are on absent contigs. The other 7 are on five contigs that are in the BAM, 17.8-63.7 kb from a contig end, so "omits 75 contigs" alone is not the cause (FC23-FC25). |

The track's evidence tables agree with each other, apart from three labels: two rows marked "agrees" carry differing values (17 vs 19 bins; 18,759 vs 17,767 bp), "occurs in no read" sits beside a 1/20 plain-join k-mer hit, and the Figure C1 caption gives n = 223 against 225 bins in the table (FC29, FC30, FC50). Limits of these re-tests: (b) reused the track's junction 31-mers; hap2 (claimed 0/16) and the unaligned contig ends were not re-tested.

## Tool status and novelty (GitHub API, 12 November 2026)

| Claim | Status |
|---|---|
| minimap2 2.31 is the latest release and uses random minimizers | ✓ v2.31 (2026-05-19); master is identical to v2.31 |
| PR lh3/minimap2#1361: plain mod-minimizer, no persisted index | ✓ open, last updated 2026-08-26; flagged as a draft only in its description, not on GitHub. Sketch-time only by design, with an explicit refusal when an index is saved or loaded. F6 attributes the silent-misuse hazard to the PR, but it came from the author's earlier prototype (FC09) |
| Qtip unmaintained | ✓ last commit 2017-11-28; "only learned predictor" ✗ (FC32) |
| DRAGEN multigenome is proprietary; 64 haplotypes; −54.2% SNV+indel errors on HG001-HG007 | ✓ (full text) |
| Table 2.1 versions: BWA-MEM2 2.3, strobealign 0.20.0, Bowtie2 2.5.5, vg 1.78.0, hifiasm 0.25.0, Verkko 2.3.2 | ✓ all current |
| HPRC MC v2.1 "-eval" graphs hold out HG002, HG005 and NA19240 | ✓ v2.1 README (HG00272 is also held out of all graphs) |
| No assembler uses methylation for phasing | ✓ for assemblers; HapBridge (bioRxiv 2025.11.07.687303) and MethPhaser already do it for reads, so the adoption gap G = 5 for U17 is overstated (FC36) |

## Citations

Parts I and III contain 40 citation links to 32 DOIs. A random sample of 40 Part II citations (seed 20261112) added 34 more DOIs, for 66 unique DOIs in all. All 66 resolve (64 via Crossref, 2 via DataCite/doi.org), with no retraction flags. All first-author surnames and years match the link text and the reference-list entries. The abstracts support the citing sentences, except Langmead 2017 ("only predictor") and Prodanov et al. 2025 (Li-Stephens). Full texts confirmed DRAGEN's 54.2% and 64 haplotypes, mapquik's Q0/Q60-only mapeval, HGSVC3's 28.3% and 1,246 centromeres, the 16.9x median of the 1,019-genome study, and the 8.4-h HG002 run for hifiasm ONT mode (in the preprint table). Two Part II citations outside the sample also fail to support their sentences: Stergachis et al. 2020 (a Fiber-seq methods paper) is cited for HPRC2 Fiber-seq data, and Li 2018 (minimap2) for string graphs, where Li 2016 (miniasm) is the right source (FC51, FC52). Part II cites both the preprint and the published versions of the Q100 benchmark and of hifiasm ONT mode (FC44).

## Consistency and cross-references

- The case studies are chapter 14 today. "Part III (chapters 11-15)", ten references to "chapter 15" and "awaits the assembly pilot (chapter 14)" are right only after the L4 chapter is inserted (FC46-FC48). All other chapter references (3, 9, 10, 11, 13) resolve correctly.
- F7 is referenced in the executive summary, Table 3.1, the roadmap and the Part III intro ("F1-F7"), but there is no F7 section, so "seven flagships" does not hold yet (FC49).
- U1-U26 match ranked_proposals.tsv exactly. Chapter 10's test labels L1-L6 and A1-A5, and chapter 9's idea labels P1-P4, collide with the pilot, source and assembly-review IDs used elsewhere (FC45).
- "232 genomes" (HPRC release-2 index) and "462 haplotypes" (QC tables; HG002 excluded) are both correct, but the report never explains the difference (FC22). "227 of them 1kGP samples" checks out against the 1kGP sample lists.

## Scope of this check

All quantitative statements in Parts I and III were traced to the track outputs (about 1,320 numeric tokens, table cells included). Design targets in the litmus criteria were checked only for consistency between chapters. The "23 further proposals" from the pilots cannot be verified until the L4 chapter exists (FC13). Spot-check scripts and outputs: [factcheck_spotchecks.tar.gz](factcheck_spotchecks.tar.gz).
