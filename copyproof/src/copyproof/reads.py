"""Read inputs for the engine: FASTQ/FASTA files (optionally gzipped) or BAM/CRAM via samtools."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from dataclasses import dataclass, field

from .engine import engine_path

_ALIGNED = (".cram", ".bam", ".sam")


@dataclass
class ReadSource:
    """One sample's reads. BAM/CRAM inputs are decoded with `samtools fastq` and streamed;
    secondary and supplementary records are dropped, unmapped reads are kept."""

    paths: list[str]
    reference: str | None = None  # FASTA that decodes the CRAMs
    samtools_threads: int = 4
    extra: list[str] = field(default_factory=list)

    def is_aligned(self) -> bool:
        kinds = {p.lower().endswith(_ALIGNED) for p in self.paths}
        if len(kinds) > 1:
            raise ValueError("mix of FASTQ/FASTA and BAM/CRAM inputs in one read set")
        return kinds.pop()

    def run(self, engine_args: list[str], log=None) -> None:
        """Run `cpk <engine_args> <reads>` over this read set."""
        cpk = engine_path()
        args = [str(a) for a in engine_args]
        if not self.is_aligned():
            subprocess.run([cpk] + args + [str(p) for p in self.paths], check=True, stderr=log)
            return
        if not shutil.which("samtools"):
            raise FileNotFoundError("samtools is needed to read BAM/CRAM inputs")
        parts = []
        for p in self.paths:
            cmd = ["samtools", "fastq", "-@", str(self.samtools_threads), "-F", "0x900"]
            if self.reference:
                cmd += ["--reference", self.reference]
            parts.append(" ".join(shlex.quote(c) for c in cmd + [str(p)]))
        producer = "( " + " ; ".join(parts) + " )"
        consumer = " ".join(shlex.quote(c) for c in [cpk] + args + ["-"])
        # pipefail: a failing samtools must fail the whole step, not truncate the counts silently
        subprocess.run(["bash", "-o", "pipefail", "-c", f"{producer} 2>/dev/null | {consumer}"], check=True, stderr=log)
