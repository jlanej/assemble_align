# How copyproof works, as built

The design in [GUIDE.md](../GUIDE.md) proposed five evidence channels around existing k-mer
tools. This is what the data led to. Where it departs from the guide, the reason is given and the
experiment is in [NOTEBOOK.md](../NOTEBOOK.md).

## 1. The k-mer engine (`cpk`)

Every channel needs the same two things: how often each k-mer of the assembly occurs in the reads,
laid out along the assembly, and which read k-mers the assembly lacks. `cpk` provides both.

**Tracked k-mers.** A canonical 31-mer is tracked when a fixed hash of it is a multiple of 4
(`--sample`). The choice depends only on the k-mer, so the same k-mers are tracked in every read
set and every assembly, and a tracked k-mer's counts are exact. Tracking a quarter of all k-mers
costs little: neighbouring k-mers are covered by the same reads and carry nearly the same
information, a heterozygous site still makes about 8 tracked k-mers, and a junction about 7.

**`cpk count`** streams a read set once. A k-mer enters the exact table at its second sighting;
the first is remembered in a Bloom filter, which keeps the singleton error k-mers (most distinct
k-mers of a read set) out of memory. Counts are exact, except that a Bloom false positive adds 1.
Reads are cut at the Illumina adapter start, because adapter read-through otherwise makes solid
k-mers that are in no genome. The database is independent of any assembly.

**`cpk scan`** loads every FASTA of a release, records each tracked k-mer's multiplicity per
*group* (haplotype 1, haplotype 2, other), joins the read databases, and writes one record per
tracked k-mer position: multiplicities and read counts. It also writes the read k-mers of the
sample that occur nowhere in the release, with their counts in the other read sets (parents).

**`cpk recruit`** and **`cpk locate`** serve the missing-sequence channel: the first keeps reads
that carry k-mers of a set, the second lists where the k-mers of some sequences occur in the
release.

Checked against a naive counter on simulated reads (exact) and against jellyfish 2.3.1 on 4 million
real reads: 13,534,076 of 13,534,369 tracked k-mers identical, 292 higher by one (Bloom false
positives), one saturated at 65,535.

## 2. Depth model

In a diploid release most tracked k-mers occur once in each haplotype. Those k-mers are almost
all correct two-copy sequence, so the median of their read counts, halved, is the **depth per
copy**. This conditions on the assembly and needs no fitting of spectrum peaks; chromosome X and Y
in a male need no special case, because their k-mers have multiplicity 1 in the release and are
expected at one copy's depth.

Depth depends on GC. The curve is the median depth of two-copy k-mers per 1% GC bin of 1-kb
windows.

**Noise.** For a window, k-mers are grouped by multiplicity m (1 to 8). The sum S of a class's
read counts has

    Var(S) = omega_m x mu x V + (psi x E)^2

- `mu` is the expected count of one k-mer and `E` the expected sum.
- `V` handles the fact that neighbouring k-mers are not independent: two k-mers d bases apart share
  reads, and their counts covary by about depth x (1 - d / R), with R the number of k-mer starts
  per read. V is the sum of (1 - d / R) over all pairs of the class's k-mers in the window closer
  than R. It is computed exactly, so the model is right for any sampling rate, read length and
  clustering of k-mers (heterozygous k-mers come in clusters of up to 31).
- `omega_m` is extra-Poisson noise, fitted per multiplicity class where the class is sparse.
- `psi` is the slow relative wander of depth that GC does not explain, fitted from the dense
  two-copy class.

Counts are winsorised per k-mer at the depth of m + 4 copies, and the noise is Student-t with 5
degrees of freedom.

## 3. S1: copy number

**Hypotheses per window.**

- *Additive, d = +1, -1, +2.* The genome holds d more copies of the window's sequence than the
  release. A k-mer of multiplicity m is expected at depth x (m + d). A copy that is missing or
  extra is not always identical to those that remain, so a change may involve every k-mer class
  or only those the release holds at least 2, 3 or 5 times (the k-mers shared between copies).
  Which classes take part is a property of the event, so it is part of the hidden state.
- *Multiplicative, factor 0, 0.25, 0.5, 0.75, 1.25, 1.5, 2.* Every class is scaled alike. This
  is uneven sequencing depth, or sequence present in a fraction of cells.

A hidden Markov model over 1-kb windows (a concordant state plus the 19 above; an event starts
with prior 10^-6 per window) gives the state path. A segment is a run of windows that *test*
their event state: where a window holds no k-mer of the classes an event involves, the event's
likelihood equals the concordant one, and such windows do not extend a segment. Segments need 20
nats over the concordant path.

Windows are not assessed when they have fewer than 30 usable k-mers, when more than half of their
k-mers exceed the multiplicity cap (8), when GC predicts less than half the depth, or when they
are mostly N; assessable islands shorter than 5 windows between such windows are dropped too.
The complement is `cannot_assess.bed`.

**Classification per segment**, from the summed evidence of its windows. Depth wanders by region
beyond what GC predicts (7.7% at 36x), but all k-mers of a stretch share that wander. So every
hypothesis has one regional factor g for all classes (prior: that width around 1, at most two
standard deviations), and each class may then differ from its expectation by 4% plus sampling
noise:

- concordant: class m at g x m;
- whole-copy change d: class m at g x (m + d) if it follows the change, g x m if it does not;
- depth factor r: class m at r x m, for any r between 0 and 3.

A hypothesis *fits* when the classes that hold 90% of the segment's k-mers, the commonest class
among them, are within 3.5 sd of it. Minor classes may not veto a fit, because they are
mixtures: a k-mer that lies twice in a surplus copy loses two copies where its classmates lose
one. They still count in the likelihood. A whole-copy change must also beat the concordant
hypothesis by 3 nats, be followed by at least half of the k-mers, and put the commonest class
within 0.35 of a whole number of copies.

| Kind | Rule | Reported as |
|---|---|---|
| `CN` | a whole-copy change fits, and is more likely than any factor that fits by 3 nats (if no factor fits: not less likely than the best one by 3 nats) | issue, high confidence |
| `CN_OR_DEPTH` | both fit: typically one k-mer class only, where one copy fewer and half the depth are the same observation | issue of low confidence if supported (below); otherwise not reported |
| `ABSENT` | every class at a quarter of its depth or less, over at least 20 windows | issue (`UNSUPPORTED`) |
| `DEPTH` | a factor fits and a whole-copy change does not | depth anomaly |
| `IRREGULAR` | neither fits | depth anomaly |

What separates a copy change from a depth change is k-mers of another multiplicity in the same
stretch. At a collapsed duplicate the k-mers unique to one haplotype sit at 2 copies where a
factor of 1.5 would put them at 1.5.

**The other copies.** An extra copy that is identical to another one has no such k-mers: the
release holds its sequence three times, the reads show two copies, and two thirds of the depth
would look the same. But its k-mers occur elsewhere in the release, and the two explanations say
different things about those places. Low sequencing depth is a property of the sequence, so it
lowers every copy, including the k-mers unique to another copy that lie between the shared ones.
Under a surplus copy those unique k-mers are read at full depth. For every `CN` and
`CN_OR_DEPTH` segment, `cpk locate --regions` lists where its k-mers occur; hits on one diagonal
form a copy; the one-copy k-mers of the windows inside those copies are tested at one copy
against the segment's own depth ratio. With at least 20 such k-mers and 3 nats either way
(their evidence alone, and together with the segment's own), an ambiguous segment becomes `CN` or
`DEPTH`. The copies are written into the issue (`other_copies`), so the sequences involved in an
event are named.

**Ambiguous segments that are reported** (low confidence), by what supports them:

- `own_kmers`: k-mers found only in the segment's own haplotype follow the change;
- `other_channel`: a missing allele, haplotype or junction, or a hapmer switch, lies within 1 kb;
- `long_and_deep`: at least 8 assessed windows, and the whole-copy reading changes the depth of
  the commonest class 1.5-fold or more. The truth genome has no such segment in 5.3 Gb
  (against 105 dips to half depth that last one or two windows), so these are rarely dips;
- `unlinked_join` (with `--long-read-bam`): the segment lies within 5 kb of a run of N that no
  long read links. A copy too many right at a scaffold join is what a false duplication at a
  join looks like, and the join in turn becomes a high-confidence issue.

**Short calls.** A whole-copy call on fewer than 5 assessed windows is of low confidence. The
typical case is 1 to 3 kb where the k-mers of one haplotype have no reads while the shared ones
sit at one copy. The truth genome has 13 of these, and an ordinary assembly of the same person
has 18 of which the truth confirms one: mostly sequence that short reads fail on in one allele.
Long reads settle them (section 6).

**Which haplotype.** The comparison is diploid: when one haplotype has collapsed a duplicate, the
shared k-mers of the homologous stretch on the other haplotype are raised too. k-mers found only
in a sequence's own group decide: if they follow the change, the issue's role is `locus`; if they
are concordant, `partner`; with too few of them, `unresolved` (for example an extra sequence and
the copy it repeats, which no k-mer can tell apart).

**Unsupported k-mers.** Runs of assembly k-mers with zero reads inside normally covered sequence
are written to `unsupported_clusters.tsv`. They are what a false join or a wrong base looks like,
and also what about a thousand correct stretches of a truth genome look like, because some
sequence is hardly read by the short-read chemistry. A run is an issue only with a second
witness: a missing allele from S2 at the same place (`BASE_ERROR`), or long reads that do not run
through it either (`FALSE_JOIN`, section 6).

## 4. S2: read sequence the release lacks

Solid read k-mers absent from the release are mostly artefacts, in four groups that had to be
removed before anything else showed:

1. adapter read-through (removed at counting);
2. k-mers spanning a homopolymer or short tandem repeat, where polymerase slippage makes reads
   with one unit more or fewer;
3. single-base errors of k-mers that occur thousands of times (satellites, rDNA, the
   mitochondrial genome, EBV): each of the 93 neighbours of such a k-mer is seen often enough to
   look solid;
4. recurrent miscalls in fixed sequence contexts (mostly T to G), which come from one read
   strand only.

Candidates are solid (at least half the depth per copy), below 3.5 copies' depth and not
low-complexity. Reads carrying at least two candidates are recruited in one pass; reads whose
k-mers are mostly high-copy in the release (by multiplicity or by read depth) are skipped, and so
are hits next to a high-copy k-mer. The other k-mers of each read are located in the release,
co-linear anchors are chained, and each read is classified by its flanks:

- both flanks close together on one sequence: `MISSING_ALLELE`, with the difference in length;
- flanks far apart or on different sequences: `MISSING_JUNCTION`, with both ends;
- one flank only: `MISSING_SEQUENCE`. This is the weakest evidence of the channel. With fewer
  than 6 reads or 5 tracked k-mers about nine in ten are artefacts where truth is known (and all
  of them on a truth genome), so those stay in the cluster table and are not issues; the others
  are low-confidence issues unless a parent's reads hold the missing k-mers too. With parents,
  one kind of weak cluster is kept: missing k-mers that are solid in the child and in neither
  parent. Recurrent artefacts recur in the parents' reads, so these are few (about 45 events
  in HG01943), and a de novo junction is among them;
- no flank: unplaced.

Placements at one locus form a cluster (clustering by shared k-mers chains unrelated loci
together through repeats). A locus is kept only if reads of both strands support it. Three or
more missing alleles less than 10 kb apart become one `MISSING_HAPLOTYPE`: one haplotype of that
stretch is not in the release. Where the flanks of an event exist at several places (both
haplotypes, or paralogous copies), it is listed at each: once as the `locus` (the place another
channel points to, if any) and elsewhere with role `twin` and a pointer `twin_of`.

A false adjacency does not make missing k-mers when the other haplotype is right: the reads' true
adjacency is still in the release. It shows as unsupported assembly k-mers (section 3).

## 5. S3: parental hapmers

A hapmer is an assembly k-mer present in the child's reads and in exactly one parent's (count at
least 3 in one, 0 in the other). Hapmers are counted in sites (runs less than k apart), because
one variant makes up to k of them. A two-state model per sequence decodes the parent of each
stretch (2% of sites assumed wrong, switch prior 10^-3 per kb). Stretches that follow the other
parent than the rest of their sequence, with at least 4 sites and 3:1, are `PHASE_SWITCH`; a
gap-delimited piece that does so is `MISPLACED_PIECE`. Sequences with fewer than 20 sites are not
assessable.

## 6. S4: long-read audit (optional)

With `--long-read-bam` (an indexed BAM of long reads against the release, local or by URL; only
slices are fetched) two things are checked.

**Joins.** For every run of at least 10 N, alignments within 10 kb are fetched. A read links the
join if one alignment runs 5 kb past it on both sides, or two alignments of the read sit on
either side in the same orientation. A join with no link and at least 3 reads within 5 kb on each
side is `UNSUPPORTED_JOIN`. That alone is ordinary, because a scaffold gap is where reads did
not connect (34 of the 42 joins of HG01943 in HPRC release 2), so the issue is of low confidence
unless another issue lies within 20 kb: a piece on the wrong haplotype, a copy too many or too
few, a missing junction. Secondary alignments are kept, because where two copies are identical
half the reads are secondary. The BAM header is compared with the FASTA first: a sequence absent
from the alignment file is "not assessable", not "unsupported".

**Assembly k-mers that no short read holds.** Two kinds of place are checked: the runs of
unsupported k-mers with normally covered flanks (the 3,000 longest), and the copy-number losses
that sit on a haplotype's own k-mers. The alignment file is used only to find the long reads at
each place. The question is then about content: are the assembly's unread k-mers in at least two
long reads? An alignment can run through a stretch the read does not hold; its k-mers cannot.

- *Carried* (at least 80% of the k-mers): the sequence is real in the long reads' DNA, and short
  reads fail on it. A run is then no issue; a short copy-number loss is dropped; a long one is
  kept at low confidence (sequence in one DNA sample and not in the other).
- *Not carried* (at most 20%, with at least 5 long reads there): the release's sequence is in no
  read of either kind. With long reads running through the place it is a `BASE_ERROR` (a wrong
  base or small consensus error); with none, a `FALSE_JOIN`. A copy-number loss stands, with
  high confidence.
- *No long read's sequence there.* At a wrong join the reads' primary alignments often go to the
  other haplotype, and secondary alignments carry no sequence. Then alignment decides: no long
  read runs through the run and at least 3 are clipped at it: `FALSE_JOIN`.

On the truth genome, 13 stretches of 1 to 3 kb have one haplotype's own k-mers unread by short
reads; HiFi reads carry all of those k-mers at 10, and carry few or none at 2, both of which are
known issues of that assembly.

## 7. Evaluation without alignment

When a truth genome of the same individual exists, its k-mers label every tracked k-mer of any
other assembly: copies in the truth minus copies in the assembly (`copyproof/truth.py`,
`scripts/eval_truth.py`). Windows where most k-mers differ are true discordances. This needs no
alignment, so the labels do not inherit an aligner's choices in repeats.
