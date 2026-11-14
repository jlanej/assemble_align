# Engine validation: are the k-mer counts right, and what do they cost?

**Question.** copyproof's channels all rest on one custom k-mer engine (`cpk`). No k-mer counter was
installed on the machine, and full counting with meryl or KMC needs hundreds of GB of temporary
space and about an hour per read set. Are a quarter of the k-mers, counted exactly in one pass,
correct and cheap enough (GUIDE decisions D8 and D9)?

## Correctness

1. **Against a naive counter** (`tests/test_engine.py`, runs on every change). A 40-kb simulated
   diploid genome with a gap, an extra sequence and 600 bp missing from both haplotypes, 20x reads
   per haplotype. Exact agreement for the read counts, for every record of the positional track
   (position, multiplicity per group, count), for the list of read k-mers absent from the release,
   for recruited reads and for anchor positions and strands.
2. **Against jellyfish 2.3.1** on the first 4,000,000 reads of the HG002 NovaSeq R1 file
   (0.60 Gb), k = 31, canonical:

   | | tracked k-mers |
   |---|---|
   | jellyfish k-mers with count ≥ 2 | 54,140,460 |
   | of those, tracked (hash multiple of 4) | 13,534,369 (24.9986%) |
   | present in the cpk database | 13,534,369 |
   | identical count | 13,534,076 |
   | cpk higher by 1 (Bloom false positive at first sighting) | 292 |
   | saturated at 65,535 (jellyfish: 73,122) | 1 |
   | in cpk only, with count 2 (singletons that hit a Bloom false positive) | 2,192 |

   No k-mer is missing and no count is lower than jellyfish's.
3. **Internal consistency.** `cpk count --exact` on hg002v1.1 and `cpk scan` of the same FASTA
   find the same 645,389,140 distinct tracked k-mers.

## Cost

Apple M1 Max, 10 cores, 64 GB; k = 31; 1 in 4 k-mers tracked.

| Step | Input | Wall time | Peak RSS |
|---|---|---|---|
| `cpk count`, uncompressed FASTQ | 20,000,000 reads, 3.02 Gb | 21 s | 19.3 GB |
| `cpk count`, one gzipped FASTQ | 355,791,257 reads, 53.72 Gb (HG002 NovaSeq R1) | 213 s | 20.4 GB |
| `jellyfish count` (all k-mers), for scale | 4,000,000 reads, 0.60 Gb | 34 s | 6.5 GB |
| `cpk count --exact` | hg002v1.1, 6.00 Gb | 27 s | 16.8 GB |
| `cpk scan`, one read set | hg002v1.1, 6.00 Gb; 1.49 billion track records | 93 s | 17.3 GB |
| window statistics, 1 process / 6 processes | 5,999,448 windows | 184 s / 71 s | |

Counting one gzipped file is bound by decompression in a single reader thread; two files are read
in parallel. The table sizes are fixed by `--genome-size` and `--sample`, so memory does not grow
with depth.

## Decisions

- **D8 (reuse or build): build.** A hash-sampled exact counter does in minutes and 20 GB what full
  counting does in an hour and hundreds of GB of scratch, and it gives positional lookups directly.
- **D9 (cost): measured.** Counting dominates only when a read set is new; a database is reused for
  every assembly of the sample.
- **D2 (k): 31** for every channel so far. One database then serves copy number, hapmers and
  junctions.

## Commands

```bash
cd copyproof/engine && cargo test --release
```

```bash
cd copyproof && .venv/bin/python -m pytest tests/test_engine.py
```

The jellyfish comparison used `jellyfish count -m 31 -C -s 1G -t 10`, `jellyfish dump -c -L 2`, and
`cpk dump`; the k-mers were encoded and hashed with numpy to select the tracked ones.
