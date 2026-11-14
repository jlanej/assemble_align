import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from copyproof import sim  # noqa: E402
from copyproof.engine import engine_path  # noqa: E402


@pytest.fixture(scope="session")
def cpk():
    return engine_path()


@pytest.fixture(scope="session")
def small_diploid(tmp_path_factory):
    """A 40-kb diploid genome and 20x reads per haplotype. The assembly differs from the genome:
    haplotype 1 has a 100-N gap, haplotype 2 ships an extra copy of 2 kb as a second sequence, and
    both haplotypes lack the 600 bp at 30,000-30,600 (sequence present only in the reads)."""
    d = tmp_path_factory.mktemp("small")
    rng = np.random.default_rng(7)
    a, b = sim.diploid(40_000, rng, het=0.002)
    reads = sim.simulate_reads({"a": a, "b": b}, depth=20, rng=rng, err=0.002)
    a_gapped = a[:15_000] + "N" * 100 + a[15_100:30_000] + a[30_600:]
    b_asm = b[:30_000] + b[30_600:]
    sim.write_fasta(d / "hap1.fa", {"h1_chr": a_gapped})
    sim.write_fasta(d / "hap2.fa.gz", {"h2_chr": b_asm, "h2_extra": b[5_000:7_000]})
    half = len(reads) // 2
    sim.write_fastq(d / "r1.fq", reads[:half])
    sim.write_fastq(d / "r2.fq.gz", reads[half:])
    return {
        "dir": d,
        "hap1": {"h1_chr": a_gapped},
        "hap2": {"h2_chr": b_asm, "h2_extra": b[5_000:7_000]},
        "reads": reads,
        "truth": {"a": a, "b": b},
    }
