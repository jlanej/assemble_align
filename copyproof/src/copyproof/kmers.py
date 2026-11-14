"""Pure-Python k-mer helpers that mirror the engine (engine/src/kmer.rs).

They are slow and exist for tests, for small sequences (junction k-mers, evidence cards) and to
document exactly which k-mers the engine tracks: a canonical k-mer is tracked when
``fmix64(kmer) % sample == 0``.
"""

from __future__ import annotations

_MASK64 = (1 << 64) - 1
_M1 = 0xFF51AFD7ED558CCD
_M2 = 0xC4CEB9FE1A85EC53
_CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def fmix64(x: int) -> int:
    """MurmurHash3 finaliser (a bijection on 64-bit integers)."""
    x ^= x >> 33
    x = (x * _M1) & _MASK64
    x ^= x >> 33
    x = (x * _M2) & _MASK64
    x ^= x >> 33
    return x


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def encode(kmer: str) -> int:
    """2-bit encoding of a k-mer (A=0, C=1, G=2, T=3), first base most significant."""
    x = 0
    for ch in kmer:
        x = (x << 2) | _CODE[ch]
    return x


def decode(x: int, k: int) -> str:
    out = []
    for _ in range(k):
        out.append("ACGT"[x & 3])
        x >>= 2
    return "".join(reversed(out))


def canonical(kmer: str) -> int:
    """Encoded canonical form: the smaller of the k-mer and its reverse complement."""
    return min(encode(kmer), encode(revcomp(kmer)))


def canonical_kmers(seq: str, k: int):
    """Yield (start, canonical code, forward_is_canonical) for each k-mer made only of A/C/G/T."""
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    fw = rv = valid = 0
    for i, ch in enumerate(seq):
        c = _CODE.get(ch)
        if c is None:
            valid = 0
            continue
        fw = ((fw << 2) | c) & mask
        rv = (rv >> 2) | ((3 - c) << shift)
        valid += 1
        if valid >= k:
            yield (i + 1 - k, fw, True) if fw < rv else (i + 1 - k, rv, False)


def tracked_kmers(seq: str, k: int, sample: int):
    """Yield (start, canonical code) for the k-mers the engine tracks at this sampling rate."""
    m = sample - 1
    for pos, code, _ in canonical_kmers(seq, k):
        if fmix64(code) & m == 0:
            yield pos, code
