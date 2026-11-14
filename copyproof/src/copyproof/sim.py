"""Small synthetic genomes and reads, for tests and for checking each evidence channel."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def random_sequence(n: int, rng: np.random.Generator, gc: float = 0.41) -> str:
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return _BASES[rng.choice(4, size=n, p=p)].tobytes().decode()


def mutate(seq: str, rate: float, rng: np.random.Generator) -> tuple[str, np.ndarray]:
    """Substitute bases at `rate`; returns the new sequence and the changed positions."""
    arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    pos = np.flatnonzero(rng.random(len(arr)) < rate)
    for p in pos:
        choices = [b for b in b"ACGT" if b != arr[p]]
        arr[p] = choices[rng.integers(3)]
    return arr.tobytes().decode(), pos


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def simulate_reads(
    seqs: dict[str, str],
    depth: float,
    rng: np.random.Generator,
    read_len: int = 150,
    err: float = 0.002,
) -> list[str]:
    """Single-end reads at `depth` per sequence copy, both strands, substitution errors at `err`."""
    reads: list[str] = []
    for seq in seqs.values():
        L = len(seq)
        if L < read_len:
            continue
        n = rng.poisson(depth * L / read_len)
        starts = rng.integers(0, L - read_len + 1, size=n)
        strands = rng.random(n) < 0.5
        for s, rc in zip(starts, strands):
            r = seq[s : s + read_len]
            if "N" in r:
                continue
            if err > 0:
                k = rng.binomial(read_len, err)
                if k:
                    b = bytearray(r.encode())
                    for p in rng.integers(0, read_len, size=k):
                        b[p] = [x for x in b"ACGT" if x != b[p]][rng.integers(3)]
                    r = b.decode()
            reads.append(revcomp(r) if rc else r)
    order = rng.permutation(len(reads))
    return [reads[i] for i in order]


def write_fasta(path, seqs: dict[str, str], width: int = 80) -> None:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i : i + width] + "\n")


def write_fastq(path, reads: list[str]) -> None:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        for i, r in enumerate(reads):
            f.write(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n")


def diploid(n: int, rng: np.random.Generator, het: float = 0.001) -> tuple[str, str]:
    """Two haplotypes of length n that differ by SNVs at rate `het`."""
    a = random_sequence(n, rng)
    b, _ = mutate(a, het, rng)
    return a, b


def read_fasta(path) -> dict[str, str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, list[str]] = {}
    name = None
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name is not None:
                out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def ensure_dir(p) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------------------------
# A synthetic trio with a child assembly that carries one planted error of each class.
# ---------------------------------------------------------------------------------------------

def make_case(rng: np.random.Generator, L: int = 600_000, het: float = 0.001) -> dict:
    """Build a truth genome, reads and an erroneous assembly.

    The child's true haplotypes are A (maternal) and B (paternal). The child's genome carries a
    30-kb tandem duplication on A, absent from both parents (a de novo event, as in HG01943).
    The assembly differs from the truth as follows (coordinates on the assembly):

    ============  ==========================================================================
    COLLAPSE      mat 100-130 kb: the genome holds this segment twice in tandem, the assembly once
    EXPANSION     mat 200-220 kb and 220-240 kb: one true copy assembled twice in tandem
    HOMOLOG_GONE  pat lacks the 15 kb homologous to mat 320-335 kb
    SWITCH        mat 420-450 kb carries the paternal haplotype's sequence
    UNSUPPORTED   pat holds 10 kb of sequence that is not in the genome
    ============  ==========================================================================
    """
    base = random_sequence(L, rng)
    A, _ = mutate(base, het / 2, rng)
    B, _ = mutate(base, het / 2, rng)
    A2, _ = mutate(base, het, rng)  # mother's other haplotype
    B2, _ = mutate(base, het, rng)  # father's other haplotype

    dup = (100_000, 130_000)
    A_truth = A[: dup[1]] + A[dup[0] : dup[1]] + A[dup[1] :]

    # assembly, maternal haplotype: start from A (which holds the duplicated segment once)
    exp = (200_000, 220_000)
    sw = (400_000, 430_000)
    mat = A[: exp[1]] + A[exp[0] : exp[1]] + A[exp[1] : sw[0]] + B[sw[0] : sw[1]] + A[sw[1] :]
    shift = exp[1] - exp[0]  # 20 kb inserted before everything downstream on mat

    # assembly, paternal haplotype
    gone = (300_000, 315_000)
    junk_at = 500_000
    junk = random_sequence(10_000, rng)
    pat = B[: gone[0]] + B[gone[1] : junk_at] + junk + B[junk_at:]
    junk_start = junk_at - (gone[1] - gone[0])

    truth = [
        ("mat", dup[0], dup[1], "COLLAPSE"),
        ("mat", exp[0], exp[1] + shift, "EXPANSION"),
        ("mat", gone[0] + shift, gone[1] + shift, "HOMOLOG_GONE"),
        ("mat", sw[0] + shift, sw[1] + shift, "SWITCH"),
        ("pat", junk_start, junk_start + len(junk), "UNSUPPORTED"),
    ]
    return {
        "child_truth": {"A": A_truth, "B": B},
        "mother": {"A": A, "A2": A2},
        "father": {"B": B, "B2": B2},
        "asm": {"mat": mat, "pat": pat},
        "truth": truth,
        "bystander": [("pat", dup[0], dup[1], "COLLAPSE"), ("pat", exp[0], exp[1], "EXPANSION")],
    }


def write_case(case: dict, outdir, rng: np.random.Generator, depth: float = 15.0, err: float = 0.002) -> dict:
    """Write the case's assembly FASTAs and the trio's FASTQs; returns the paths."""
    d = ensure_dir(outdir)
    paths = {"mat": d / "asm.mat.fa", "pat": d / "asm.pat.fa"}
    write_fasta(paths["mat"], {"mat": case["asm"]["mat"]})
    write_fasta(paths["pat"], {"pat": case["asm"]["pat"]})
    for who, key in (("child", "child_truth"), ("mother", "mother"), ("father", "father")):
        reads = simulate_reads(case[key], depth=depth, rng=rng, err=err)
        paths[who] = d / f"{who}.fq.gz"
        write_fastq(paths[who], reads)
    with open(d / "truth.bed", "w") as f:
        for c, s, e, label in case["truth"]:
            f.write(f"{c}\t{s}\t{e}\t{label}\n")
    paths["truth"] = d / "truth.bed"
    return paths
