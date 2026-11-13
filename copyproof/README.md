# copyproof

*Proofreading released genome assemblies for copies that are missing, doubled, misjoined or on the wrong haplotype.*

**Status: design.** This directory holds the plan for the tool. No code has been written yet.

| File | What it is |
|---|---|
| [GUIDE.md](GUIDE.md) | What to build, why, with which data, how to tell whether it works, and the decisions to take along the way |
| [NOTEBOOK.md](NOTEBOOK.md) | Running log of experiments and decisions, with the baselines to reproduce |
| [resources.tsv](resources.tsv) | Public data and in-repository assets the work uses (39 entries, checked 2026-11-13) |

## The idea

Assembly QC in current releases asks whether read coverage looks wrong. Some copies are real but misplaced: scaffolded into the wrong haplotype, collapsed onto a twin copy, or duplicated at a join. Their coverage is normal, so they pass.

copyproof asks other questions of every released haplotype:

- Does each sequence occur as many times in the assembly as the sample's reads say it should?
- Are any read k-mers missing from the assembly, and do they mark a missing junction?
- Does every scaffold piece carry its own parent's sequence?
- Does any read span each join?

It combines the answers into one evidence card per issue and marks what it cannot assess. It is meant to run across whole releases (HPRC release 2, HGSVC3) and to produce a curated issue list and an estimate of how common the problem is.

The motivating cases, HG01943 and HG01103, are in [`../litmus/cases_HG01943_HG01103/`](../litmus/cases_HG01943_HG01103/). The proposal is F4 in the report ([`../report/assemble_align_review.md`](../report/assemble_align_review.md)).

## Picking this up

1. Read GUIDE §1–§3 (the idea and the evidence), §8 (the validation ladder) and §10 (milestones).
2. Milestone M0: re-run the two case drivers and confirm the baselines listed in NOTEBOOK.md.
3. Then M1 (the copy-number channel) and M2 (the early false-flag test on HG002), which decides whether the approach scales.
4. Log every experiment in NOTEBOOK.md: the question, the run, the numbers and the decision.

The guide is a map rather than a specification. Where results disagree with it, follow the results and record why.

## Conventions

- Reads, alignments, k-mer databases and assemblies stay out of git; the repository `.gitignore` excludes them. Scratch space comes from `mktemp` and is removed on exit.
- Flagged loci are "issues" until adjudicated. Share them with the releasing consortia as issues.
- Licences follow the repository root: code under MIT, text under CC BY 4.0.
