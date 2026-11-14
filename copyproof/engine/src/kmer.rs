//! 2-bit k-mer encoding, canonical rolling k-mers and the invertible hash used for sampling.
//!
//! A k-mer is tracked when `fmix64(canonical) & (sample - 1) == 0`, so the decision depends only on
//! the k-mer itself: a tracked k-mer is tracked in every read set and every assembly.

/// Marks an empty table slot. `fmix64` is a bijection on u64, so at most one k-mer hashes here;
/// it is treated as never tracked.
pub const EMPTY: u64 = u64::MAX;

/// A=0, C=1, G=2, T=3, anything else 4.
pub const ENC: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

const M1: u64 = 0xff51_afd7_ed55_8ccd;
const M2: u64 = 0xc4ce_b9fe_1a85_ec53;

/// Modular inverse of an odd u64 (Newton iteration).
const fn inv64(x: u64) -> u64 {
    let mut y = x;
    let mut i = 0;
    while i < 6 {
        y = y.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(y)));
        i += 1;
    }
    y
}
const M1_INV: u64 = inv64(M1);
const M2_INV: u64 = inv64(M2);

/// MurmurHash3 finaliser: a bijection on u64 with good avalanche.
#[inline(always)]
pub fn fmix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(M1);
    x ^= x >> 33;
    x = x.wrapping_mul(M2);
    x ^= x >> 33;
    x
}

/// Inverse of `fmix64`.
#[inline]
pub fn fmix64_inv(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(M2_INV);
    x ^= x >> 33;
    x = x.wrapping_mul(M1_INV);
    x ^= x >> 33;
    x
}

/// Call `f(start, canonical, forward_is_canonical)` for every k-mer of `seq` made only of A/C/G/T.
/// `seq` holds ASCII bases; any other byte (N, newline, ...) breaks the run.
#[inline(always)]
pub fn for_each_canon<F: FnMut(usize, u64, bool)>(seq: &[u8], k: usize, mut f: F) {
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = ENC[b as usize] as u64;
        if c < 4 {
            fw = ((fw << 2) | c) & mask;
            rv = (rv >> 2) | ((3 - c) << shift);
            valid += 1;
            if valid >= k {
                if fw < rv {
                    f(i + 1 - k, fw, true)
                } else {
                    f(i + 1 - k, rv, false)
                }
            }
        } else {
            valid = 0;
        }
    }
}

/// Same as `for_each_canon` for a slice of codes (0..=3, 4 = not a base).
#[inline(always)]
pub fn for_each_canon_codes<F: FnMut(usize, u64, bool)>(codes: &[u8], k: usize, mut f: F) {
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for (i, &c8) in codes.iter().enumerate() {
        let c = c8 as u64;
        if c < 4 {
            fw = ((fw << 2) | c) & mask;
            rv = (rv >> 2) | ((3 - c) << shift);
            valid += 1;
            if valid >= k {
                if fw < rv {
                    f(i + 1 - k, fw, true)
                } else {
                    f(i + 1 - k, rv, false)
                }
            }
        } else {
            valid = 0;
        }
    }
}

/// Decode a 2-bit k-mer to ASCII.
pub fn decode(mut x: u64, k: usize) -> String {
    let mut s = vec![b'A'; k];
    for i in (0..k).rev() {
        s[i] = b"ACGT"[(x & 3) as usize];
        x >>= 2;
    }
    String::from_utf8(s).unwrap()
}

/// Encode an ASCII k-mer to its canonical 2-bit form; None if it holds a non-ACGT base.
#[cfg(test)]
pub fn encode_canon(s: &[u8]) -> Option<u64> {
    let k = s.len();
    let mut out = None;
    for_each_canon(s, k, |_, c, _| out = Some(c));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_is_invertible() {
        let mut x = 0x1234_5678_9abc_def0u64;
        for _ in 0..10_000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            assert_eq!(fmix64_inv(fmix64(x)), x);
            assert_eq!(fmix64(fmix64_inv(x)), x);
        }
        assert_eq!(M1.wrapping_mul(M1_INV), 1);
        assert_eq!(M2.wrapping_mul(M2_INV), 1);
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    }

    #[test]
    fn canonical_is_strand_symmetric() {
        let seq = b"ACGTTGCAAGGCTTAACGNACGATCGATCGGGATATATCCGGAATTCCGGAAT";
        let k = 11;
        let mut a = vec![];
        for_each_canon(seq, k, |p, c, _| a.push((p, c)));
        let rc = revcomp(seq);
        let mut b = vec![];
        for_each_canon(&rc, k, |p, c, _| b.push((seq.len() - k - p, c)));
        b.sort();
        assert_eq!(a, b);
        // the N at index 18 removes every k-mer that overlaps it
        assert!(a.iter().all(|&(p, _)| p + k <= 18 || p > 18));
        assert_eq!(a.len(), (18 - k + 1) + (seq.len() - 19 - k + 1));
    }

    #[test]
    fn encode_decode_roundtrip() {
        let s = b"ACGTACGTTTGACCA";
        let c = encode_canon(s).unwrap();
        let d = decode(c, s.len());
        assert!(d.as_bytes() == s || d.as_bytes() == revcomp(s).as_slice());
        assert!(encode_canon(b"ACGNA").is_none());
    }
}
