"""Read the engine's scan outputs: the positional k-mer track and its side files."""

from __future__ import annotations

import json
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


def _read_header(path: Path, magic: bytes) -> tuple[dict, int]:
    with open(path, "rb") as f:
        got = f.read(8)
        if got != magic:
            raise ValueError(f"{path}: expected {magic!r}, found {got!r}")
        (n,) = struct.unpack("<I", f.read(4))
        header = json.loads(f.read(n).decode())
    return header, 12 + n


@dataclass
class Scan:
    """One `cpk scan` output set, addressed by its prefix."""

    prefix: str
    meta: dict = field(init=False)
    header: dict = field(init=False)
    contigs: pd.DataFrame = field(init=False)
    gaps: pd.DataFrame = field(init=False)

    def __post_init__(self):
        p = self.prefix
        self.meta = json.loads(Path(f"{p}.scan.json").read_text())
        self.header, off = _read_header(Path(f"{p}.track.bin"), b"CPKTRK01")
        self.k = int(self.header["k"])
        self.sample = int(self.header["sample"])
        self.groups: list[str] = list(self.header["groups"])
        self.reads: list[str] = list(self.header["reads"])
        self.bin = int(self.header["bin"])
        G, R = len(self.groups), len(self.reads)
        fields = [("gap", "u1"), ("m", "u1", (G,))]
        if R:
            fields.append(("c", "<u2", (R,)))
        self.dtype = np.dtype(fields)
        assert self.dtype.itemsize == self.header["record_bytes"]
        self.track = np.memmap(f"{p}.track.bin", dtype=self.dtype, mode="r", offset=off)
        self.contigs = pd.read_csv(f"{p}.contigs.tsv", sep="\t", dtype={"name": str, "group": str})
        self.gc = np.memmap(f"{p}.gc{self.bin}.u8", dtype=np.uint8, mode="r")
        self.acgt = np.memmap(f"{p}.acgt{self.bin}.u8", dtype=np.uint8, mode="r")
        gp = Path(f"{p}.gaps.bed")
        if gp.stat().st_size:
            self.gaps = pd.read_csv(gp, sep="\t", names=["contig", "start", "end"], dtype={"contig": str})
        else:
            self.gaps = pd.DataFrame({"contig": pd.Series(dtype=str), "start": pd.Series(dtype=int), "end": pd.Series(dtype=int)})
        self._by_name = {n: i for i, n in enumerate(self.contigs["name"])}

    # ---- per-sequence access -------------------------------------------------------------

    def contig_index(self, contig) -> int:
        return int(contig) if not isinstance(contig, str) else self._by_name[contig]

    def group_index(self, contig) -> int:
        return self.groups.index(self.contigs.at[self.contig_index(contig), "group"])

    def records(self, contig):
        """Tracked k-mers of one sequence: positions (int64), m (n x groups), c (n x read sets)."""
        row = self.contigs.iloc[self.contig_index(contig)]
        rec = self.track[int(row.rec_start) : int(row.rec_start) + int(row.n_rec)]
        pos = np.cumsum(rec["gap"], dtype=np.int64)
        m = np.asarray(rec["m"])
        real = m.any(axis=1)
        c = np.asarray(rec["c"])[real] if self.reads else np.zeros((int(real.sum()), 0), dtype=np.uint16)
        return pos[real], m[real], c

    def region(self, contig, start: int, end: int):
        """Tracked k-mers whose start lies in [start, end)."""
        pos, m, c = self.records(contig)
        a, b = np.searchsorted(pos, [start, end])
        return pos[a:b], m[a:b], c[a:b]

    def composition(self, contig):
        """Per-bin counts of G+C bases and of A/C/G/T bases for one sequence."""
        row = self.contigs.iloc[self.contig_index(contig)]
        a, n = int(row.bin_start), int(row.n_bins)
        return np.asarray(self.gc[a : a + n]), np.asarray(self.acgt[a : a + n])

    # ---- read k-mers absent from the release ----------------------------------------------

    def novel(self):
        """Read k-mers (first read set) absent from the release: (kmers u64, counts n x read sets)."""
        path = Path(f"{self.prefix}.novel.bin")
        header, off = _read_header(path, b"CPKNOV01")
        R = len(header["reads"])
        dt = np.dtype([("kmer", "<u8"), ("c", "<u2", (R,))])
        arr = np.memmap(path, dtype=dt, mode="r", offset=off)
        return np.asarray(arr["kmer"]), np.asarray(arr["c"])

    def spectrum(self) -> pd.DataFrame:
        """Distinct tracked k-mers by multiplicity per group (each capped at 3) and read count."""
        return pd.read_csv(f"{self.prefix}.spectrum.tsv", sep="\t")
