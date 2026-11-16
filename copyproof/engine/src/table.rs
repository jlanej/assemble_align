//! Lock-free open-addressing table of k-mer hashes, a blocked Bloom filter and a HyperLogLog.

use crate::kmer::{fmix64, EMPTY};
use std::sync::atomic::{AtomicU16, AtomicU32, AtomicU64, Ordering::Relaxed};

pub fn atomic_u64_vec(n: usize, v: u64) -> Vec<AtomicU64> {
    let raw = vec![v; n];
    // AtomicU64 has the same size, alignment and bit validity as u64.
    unsafe { std::mem::transmute::<Vec<u64>, Vec<AtomicU64>>(raw) }
}
pub fn atomic_u32_vec(n: usize) -> Vec<AtomicU32> {
    let raw = vec![0u32; n];
    unsafe { std::mem::transmute::<Vec<u32>, Vec<AtomicU32>>(raw) }
}
pub fn atomic_u16_vec(n: usize) -> Vec<AtomicU16> {
    let raw = vec![0u16; n];
    unsafe { std::mem::transmute::<Vec<u16>, Vec<AtomicU16>>(raw) }
}

/// Add `n` to a u16 counter, saturating at 65,535.
#[inline]
pub fn sat_add_u16(c: &AtomicU16, n: u32) {
    let n = n.min(65_535) as u16;
    let mut cur = c.load(Relaxed);
    loop {
        let new = cur.saturating_add(n);
        if new == cur {
            return;
        }
        match c.compare_exchange_weak(cur, new, Relaxed, Relaxed) {
            Ok(_) => return,
            Err(x) => cur = x,
        }
    }
}

/// Keys are k-mer hashes; a key's home slot comes from the hash's high bits, so the low bits
/// (used for sampling) do not bias placement. Linear probing; keys are never removed.
pub struct KeyTable {
    pub keys: Vec<AtomicU64>,
    n: u64,
}

impl KeyTable {
    pub fn new(nslots: usize) -> KeyTable {
        let n = nslots.max(16);
        KeyTable { keys: atomic_u64_vec(n, EMPTY), n: n as u64 }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.n as usize
    }

    #[inline(always)]
    pub fn home(&self, h: u64) -> usize {
        ((h as u128 * self.n as u128) >> 64) as usize
    }

    /// Slot of `h`, or None if it is absent.
    #[inline]
    pub fn find(&self, h: u64) -> Option<usize> {
        self.find_from(h, self.home(h)).ok()
    }

    /// Probe from slot `i`: Ok(slot) if found, Err(first empty slot) if absent.
    #[inline]
    pub fn find_from(&self, h: u64, mut i: usize) -> Result<usize, usize> {
        let n = self.n as usize;
        let mut probes = 0usize;
        loop {
            let k = self.keys[i].load(Relaxed);
            if k == h {
                return Ok(i);
            }
            if k == EMPTY {
                return Err(i);
            }
            i += 1;
            if i == n {
                i = 0;
            }
            probes += 1;
            if probes > n {
                panic!("k-mer table is full");
            }
        }
    }

    /// Insert `h` if absent, probing from slot `i`. Returns (slot, newly inserted).
    #[inline]
    pub fn insert_from(&self, h: u64, mut i: usize) -> (usize, bool) {
        let n = self.n as usize;
        let mut probes = 0usize;
        loop {
            let k = self.keys[i].load(Relaxed);
            if k == h {
                return (i, false);
            }
            if k == EMPTY {
                match self.keys[i].compare_exchange(EMPTY, h, Relaxed, Relaxed) {
                    Ok(_) => return (i, true),
                    Err(cur) => {
                        if cur == h {
                            return (i, false);
                        }
                    }
                }
            }
            i += 1;
            if i == n {
                i = 0;
            }
            probes += 1;
            if probes > n {
                panic!("k-mer table is full; raise --genome-size or --sample");
            }
        }
    }

    #[inline]
    pub fn insert(&self, h: u64) -> (usize, bool) {
        self.insert_from(h, self.home(h))
    }
}

/// Blocked Bloom filter: two bits inside one 512-bit block, so a query touches one cache line.
pub struct Bloom {
    words: Vec<AtomicU64>,
    nblocks: u64,
}

impl Bloom {
    pub fn new(bits: u64) -> Bloom {
        let nblocks = (bits / 512).max(1);
        Bloom { words: atomic_u64_vec((nblocks * 8) as usize, 0), nblocks }
    }

    pub fn bits(&self) -> u64 {
        self.nblocks * 512
    }

    /// Set the key's bits; true if they were all set already (the key was probably seen before).
    #[inline]
    pub fn test_and_set(&self, h: u64) -> bool {
        let g = fmix64(h ^ 0x5851_f42d_4c95_7f2d);
        let block = ((g as u128 * self.nblocks as u128) >> 64) as usize * 8;
        let b1 = (g & 511) as usize;
        let b2 = ((g >> 9) & 511) as usize;
        let (w1, m1) = (block + (b1 >> 6), 1u64 << (b1 & 63));
        let (w2, m2) = (block + (b2 >> 6), 1u64 << (b2 & 63));
        let o1 = self.words[w1].load(Relaxed);
        let o2 = self.words[w2].load(Relaxed);
        if o1 & m1 != 0 && o2 & m2 != 0 {
            return true;
        }
        let p1 = self.words[w1].fetch_or(m1, Relaxed);
        let p2 = self.words[w2].fetch_or(m2, Relaxed);
        p1 & m1 != 0 && p2 & m2 != 0
    }
}

/// HyperLogLog with 2^16 registers (about 0.4% standard error); mergeable.
pub struct Hll {
    reg: Vec<u8>,
}

const HLL_P: u32 = 16;

impl Hll {
    pub fn new() -> Hll {
        Hll { reg: vec![0u8; 1 << HLL_P] }
    }

    #[inline]
    pub fn add(&mut self, h: u64) {
        let g = fmix64(h ^ 0x2545_f491_4f6c_dd1d);
        let idx = (g >> (64 - HLL_P)) as usize;
        let rest = (g << HLL_P) | (1u64 << (HLL_P - 1));
        let rank = rest.leading_zeros() as u8 + 1;
        if rank > self.reg[idx] {
            self.reg[idx] = rank;
        }
    }

    pub fn merge(&mut self, other: &Hll) {
        for (a, b) in self.reg.iter_mut().zip(other.reg.iter()) {
            if *b > *a {
                *a = *b;
            }
        }
    }

    pub fn estimate(&self) -> f64 {
        let m = self.reg.len() as f64;
        let alpha = 0.7213 / (1.0 + 1.079 / m);
        let sum: f64 = self.reg.iter().map(|&r| 2f64.powi(-(r as i32))).sum();
        let raw = alpha * m * m / sum;
        let zeros = self.reg.iter().filter(|&&r| r == 0).count() as f64;
        if raw <= 2.5 * m && zeros > 0.0 {
            m * (m / zeros).ln()
        } else {
            raw
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_insert_find() {
        let t = KeyTable::new(1000);
        let keys: Vec<u64> = (1..=600u64).map(|i| fmix64(i * 7919)).collect();
        for &k in &keys {
            let (_, new) = t.insert(k);
            assert!(new);
        }
        for &k in &keys {
            let (i, new) = t.insert(k);
            assert!(!new);
            assert_eq!(t.find(k), Some(i));
        }
        assert_eq!(t.find(fmix64(123_456_789)), None);
    }

    #[test]
    fn table_concurrent_inserts_are_unique() {
        use std::sync::Arc;
        let t = Arc::new(KeyTable::new(200_000));
        let counts = Arc::new(atomic_u16_vec(200_000));
        let mut hs = vec![];
        for th in 0..8 {
            let (t, counts) = (t.clone(), counts.clone());
            hs.push(std::thread::spawn(move || {
                let mut newly = 0u64;
                for i in 0..100_000u64 {
                    // every thread inserts the same 100k keys, in a different order
                    let key = fmix64(((i * 31 + th * 977) % 100_000) + 1);
                    let (slot, new) = t.insert(key);
                    newly += new as u64;
                    sat_add_u16(&counts[slot], 1);
                }
                newly
            }));
        }
        let newly: u64 = hs.into_iter().map(|h| h.join().unwrap()).sum();
        assert_eq!(newly, 100_000);
        let filled = t.keys.iter().filter(|k| k.load(Relaxed) != EMPTY).count();
        assert_eq!(filled, 100_000);
        let total: u64 = counts.iter().map(|c| c.load(Relaxed) as u64).sum();
        assert_eq!(total, 800_000);
        for i in 0..100_000u64 {
            let s = t.find(fmix64(i + 1)).unwrap();
            assert_eq!(counts[s].load(Relaxed), 8);
        }
    }

    #[test]
    fn saturating_counter() {
        let c = AtomicU16::new(65_530);
        sat_add_u16(&c, 3);
        assert_eq!(c.load(Relaxed), 65_533);
        sat_add_u16(&c, 100);
        assert_eq!(c.load(Relaxed), 65_535);
        sat_add_u16(&c, 1);
        assert_eq!(c.load(Relaxed), 65_535);
    }

    #[test]
    fn bloom_remembers() {
        let b = Bloom::new(1 << 20);
        let mut fp = 0;
        for i in 0..20_000u64 {
            if b.test_and_set(fmix64(i)) {
                fp += 1;
            }
        }
        assert!(fp < 400, "false positives while filling: {fp}");
        for i in 0..20_000u64 {
            assert!(b.test_and_set(fmix64(i)));
        }
    }

    #[test]
    fn hll_estimates() {
        let mut h = Hll::new();
        for i in 0..2_000_000u64 {
            h.add(fmix64(i));
            h.add(fmix64(i)); // duplicates do not count
        }
        let e = h.estimate();
        assert!((e / 2_000_000.0 - 1.0).abs() < 0.02, "estimate {e}");
    }
}
