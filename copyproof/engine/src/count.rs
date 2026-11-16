//! `cpk count`: one streaming pass over a read set, producing a database of tracked k-mer counts.
//!
//! Only k-mers whose hash passes the sampling test are tracked. A k-mer enters the exact table at
//! its second sighting; the first sighting is remembered in a Bloom filter, which keeps the
//! singleton sequencing-error k-mers (most distinct k-mers of a read set) out of the table.
//! Counts are therefore exact, except that a Bloom false positive adds one to a k-mer's count.

use crate::fastx::{stream_reads, Batch, ReadStats};
use crate::kmer::{decode, fmix64, fmix64_inv, for_each_canon, EMPTY};
use crate::table::{atomic_u16_vec, sat_add_u16, Bloom, KeyTable};
use crate::util::{json_str, peak_rss_bytes, threads_from, Args, Res};
use crossbeam_channel::{bounded, Receiver};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, AtomicU16, AtomicU64, AtomicUsize, Ordering::Relaxed};

pub const DB_MAGIC: &[u8; 8] = b"CPKDB001";
pub const DB_ENTRY: usize = 10;
const FLUSH: usize = 64;

/// Start of the Illumina TruSeq adapters (read 1 and read 2 share it). A read from a fragment
/// shorter than the read runs into the adapter; left in, it makes k-mers that are in no genome.
pub const TRUSEQ: &[u8] = b"AGATCGGAAGAGC";

/// Cut a read at the first full occurrence of an adapter start, or at a partial occurrence of
/// at least `MIN_TAIL` bases at the read's end.
#[inline]
pub fn trim_adapter<'a>(seq: &'a [u8], finders: &[(memchr::memmem::Finder<'static>, Vec<u8>)]) -> &'a [u8] {
    const MIN_TAIL: usize = 5;
    let mut end = seq.len();
    for (finder, adapter) in finders {
        if let Some(i) = finder.find(&seq[..end]) {
            end = i;
            continue;
        }
        let max = adapter.len().min(end);
        let mut j = max.saturating_sub(1);
        while j >= MIN_TAIL {
            if seq[..end].ends_with(&adapter[..j]) {
                end -= j;
                break;
            }
            j -= 1;
        }
    }
    &seq[..end]
}

pub fn adapter_finders(a: &Args, exact: bool) -> Vec<(memchr::memmem::Finder<'static>, Vec<u8>)> {
    let given = a.all("adapter");
    let list: Vec<Vec<u8>> = if exact || given.iter().any(|x| x == "none") {
        vec![]
    } else if given.is_empty() {
        vec![TRUSEQ.to_vec()]
    } else {
        given.iter().map(|x| x.to_ascii_uppercase().into_bytes()).collect()
    };
    list.into_iter().map(|v| (memchr::memmem::Finder::new(&v).into_owned(), v)).collect()
}

pub fn check_k_sample(k: usize, s: u64) -> Res<()> {
    if !(11..=31).contains(&k) || k % 2 == 0 {
        return Err(format!("--k must be odd and between 11 and 31 (got {k})").into());
    }
    if s == 0 || !s.is_power_of_two() {
        return Err(format!("--sample must be a power of two (got {s})").into());
    }
    Ok(())
}

#[derive(Default)]
struct WorkerStats {
    kmers: u64,
    tracked: u64,
}

struct Ctx<'a> {
    table: &'a KeyTable,
    counts: &'a [AtomicU16],
    bloom: &'a Bloom,
    k: usize,
    smask: u64,
    inserted: &'a AtomicU64,
    overflow: &'a AtomicBool,
    max_keys: u64,
    /// Count every sighting (no Bloom filter): for assemblies, where count 1 is meaningful.
    exact: bool,
    adapters: &'a [(memchr::memmem::Finder<'static>, Vec<u8>)],
    trimmed: &'a AtomicU64,
}

#[inline]
fn flush(buf: &mut Vec<(u64, u32)>, ctx: &Ctx, inserted: &mut u64) {
    let n = buf.len();
    let mut home = [0usize; FLUSH];
    let mut k0 = [0u64; FLUSH];
    for j in 0..n {
        home[j] = ctx.table.home(buf[j].0);
    }
    // Independent loads first: the memory accesses of a whole batch overlap in time.
    for j in 0..n {
        k0[j] = ctx.table.keys[home[j]].load(Relaxed);
        let _ = ctx.counts[home[j]].load(Relaxed);
    }
    for j in 0..n {
        let (h, c) = buf[j];
        if k0[j] == h {
            sat_add_u16(&ctx.counts[home[j]], c);
            continue;
        }
        match ctx.table.find_from(h, home[j]) {
            Ok(i) => sat_add_u16(&ctx.counts[i], c),
            Err(empty) if ctx.exact => {
                let (i, new) = ctx.table.insert_from(h, empty);
                sat_add_u16(&ctx.counts[i], c);
                *inserted += new as u64;
            }
            Err(empty) => {
                let seen = ctx.bloom.test_and_set(h);
                if seen || c >= 2 {
                    let (i, new) = ctx.table.insert_from(h, empty);
                    // `seen` stands for one earlier sighting that only the Bloom filter recorded.
                    let add = if new && seen { c + 1 } else { c };
                    sat_add_u16(&ctx.counts[i], add);
                    *inserted += new as u64;
                }
            }
        }
    }
    buf.clear();
}

fn worker(rx: Receiver<Batch>, ctx: &Ctx) -> WorkerStats {
    let mut st = WorkerStats::default();
    let mut buf: Vec<(u64, u32)> = Vec::with_capacity(FLUSH);
    let mut inserted = 0u64;
    for batch in rx.iter() {
        if ctx.overflow.load(Relaxed) {
            continue; // drain the channel so the readers can finish
        }
        let (mut kmers, mut tracked, mut trimmed) = (0u64, 0u64, 0u64);
        let mut each = |seq: &[u8], buf: &mut Vec<(u64, u32)>, inserted: &mut u64| {
            for_each_canon(seq, ctx.k, |_, canon, _| {
                kmers += 1;
                let h = fmix64(canon);
                if h & ctx.smask == 0 && h != EMPTY {
                    tracked += 1;
                    if let Some(last) = buf.last_mut() {
                        if last.0 == h {
                            last.1 += 1; // runs of one k-mer (homopolymers) are added in one step
                            return;
                        }
                    }
                    buf.push((h, 1));
                    if buf.len() == FLUSH {
                        flush(buf, ctx, inserted);
                    }
                }
            });
        };
        if ctx.adapters.is_empty() {
            // any byte that is not a base ends a k-mer, so the batch can be scanned as one sequence
            each(&batch.data, &mut buf, &mut inserted);
        } else {
            for seq in batch.data.split(|&b| b == b'\n') {
                let cut = trim_adapter(seq, ctx.adapters);
                trimmed += (cut.len() < seq.len()) as u64;
                each(cut, &mut buf, &mut inserted);
                // a read's last k-mer must not merge with the next read's first
                flush(&mut buf, ctx, &mut inserted);
            }
        }
        flush(&mut buf, ctx, &mut inserted);
        ctx.trimmed.fetch_add(trimmed, Relaxed);
        st.kmers += kmers;
        st.tracked += tracked;
        if inserted > 0 {
            let total = ctx.inserted.fetch_add(inserted, Relaxed) + inserted;
            inserted = 0;
            if total > ctx.max_keys {
                ctx.overflow.store(true, Relaxed);
            }
        }
    }
    st
}

pub fn run(a: &Args) -> Res<()> {
    let k = a.usize_or("k", 31)?;
    let s = a.usize_or("sample", 4)? as u64;
    check_k_sample(k, s)?;
    let threads = threads_from(a)?;
    let exact = a.get("exact").map(|v| v != "false" && v != "0").unwrap_or(false);
    let g = a.usize_or("genome-size", 3_100_000_000)? as f64;
    let slots = a.usize_or("table-slots", (2.0 * g / s as f64).ceil() as usize)?.max(1 << 16);
    let bloom_bits = if exact { 1 << 16 } else { a.usize_or("bloom-bits", (36.0 * g / s as f64).ceil() as usize)?.max(1 << 16) as u64 };
    let min_count = a.usize_or("min-count", if exact { 1 } else { 2 })?.clamp(1, 65_535) as u16;
    let out = a.req("out")?.to_string();
    let inputs = a.pos.clone();
    if inputs.is_empty() {
        return Err("count: give at least one FASTQ/FASTA input ('-' for stdin)".into());
    }
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

    log!(
        "count: k={k} sample=1/{s} threads={threads} table={:.2} GB ({slots} slots) bloom={:.2} GB{}",
        slots as f64 * 10.0 / 1e9,
        bloom_bits as f64 / 8e9,
        if exact { " (exact: every sighting counted)" } else { "" }
    );
    let table = KeyTable::new(slots);
    let counts = atomic_u16_vec(slots);
    let bloom = Bloom::new(bloom_bits);
    let inserted = AtomicU64::new(0);
    let overflow = AtomicBool::new(false);
    let adapters = adapter_finders(a, exact);
    let trimmed = AtomicU64::new(0);
    let ctx = Ctx {
        table: &table,
        counts: &counts,
        bloom: &bloom,
        k,
        smask: s - 1,
        inserted: &inserted,
        overflow: &overflow,
        max_keys: (slots as f64 * 0.9) as u64,
        exact,
        adapters: &adapters,
        trimmed: &trimmed,
    };

    let nreaders = inputs.len().min((threads / 2).max(1));
    let next_input = AtomicUsize::new(0);
    let (tx, rx) = bounded::<Batch>(threads * 4);
    let mut read_stats = ReadStats::default();
    let mut wstats = WorkerStats::default();
    let mut reader_err: Option<String> = None;
    std::thread::scope(|sc| {
        let readers: Vec<_> = (0..nreaders)
            .map(|_| {
                let tx = tx.clone();
                let (inputs, next_input) = (&inputs, &next_input);
                sc.spawn(move || -> Result<ReadStats, String> {
                    let mut st = ReadStats::default();
                    loop {
                        let i = next_input.fetch_add(1, Relaxed);
                        if i >= inputs.len() {
                            break;
                        }
                        let r = stream_reads(&inputs[i], i, 4 << 20, k - 1, &tx).map_err(|e| e.to_string())?;
                        log!("count: finished {} ({} records, {:.2} Gb)", inputs[i], r.records, r.bases as f64 / 1e9);
                        st.records += r.records;
                        st.bases += r.bases;
                    }
                    Ok(st)
                })
            })
            .collect();
        drop(tx);
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let rx = rx.clone();
                let ctx = &ctx;
                sc.spawn(move || worker(rx, ctx))
            })
            .collect();
        drop(rx);
        for r in readers {
            match r.join().unwrap() {
                Ok(st) => {
                    read_stats.records += st.records;
                    read_stats.bases += st.bases;
                }
                Err(e) => reader_err = Some(e),
            }
        }
        for w in workers {
            let st = w.join().unwrap();
            wstats.kmers += st.kmers;
            wstats.tracked += st.tracked;
        }
    });
    if let Some(e) = reader_err {
        return Err(e.into());
    }
    if overflow.load(Relaxed) {
        return Err(format!(
            "the k-mer table filled beyond 90% ({} keys in {slots} slots); rerun with a larger --genome-size or --sample",
            inserted.load(Relaxed)
        )
        .into());
    }
    let nkeys = inserted.load(Relaxed);
    log!(
        "count: {} reads ({} cut at an adapter), {:.2} Gb, {:.3e} k-mers, {:.3e} tracked; {} keys (load {:.2})",
        read_stats.records,
        trimmed.load(Relaxed),
        read_stats.bases as f64 / 1e9,
        wstats.kmers as f64,
        wstats.tracked as f64,
        nkeys,
        nkeys as f64 / slots as f64
    );

    // Histogram and entry count.
    const CH: usize = 1 << 22;
    let nchunks = slots.div_ceil(CH);
    let (hist, n_entries) = (0..nchunks)
        .into_par_iter()
        .map(|ci| {
            let mut hist = vec![0u64; 65_536];
            let mut n = 0u64;
            for i in ci * CH..((ci + 1) * CH).min(slots) {
                if table.keys[i].load(Relaxed) != EMPTY {
                    let c = counts[i].load(Relaxed);
                    hist[c as usize] += 1;
                    n += (c >= min_count) as u64;
                }
            }
            (hist, n)
        })
        .reduce(
            || (vec![0u64; 65_536], 0u64),
            |mut a, b| {
                for (x, y) in a.0.iter_mut().zip(b.0.iter()) {
                    *x += *y;
                }
                (a.0, a.1 + b.1)
            },
        );

    let header = format!(
        "{{\"format\": \"cpkdb\", \"version\": 1, \"k\": {k}, \"sample\": {s}, \"n_entries\": {n_entries}, \"min_count\": {min_count}, \"exact\": {exact}, \"adapters\": [{}], \"adapter_trimmed_reads\": {}, \"records\": {}, \"bases\": {}, \"kmers\": {}, \"tracked\": {}, \"table_slots\": {slots}, \"table_keys\": {nkeys}, \"bloom_bits\": {}, \"inputs\": [{}]}}",
        adapters.iter().map(|(_, v)| json_str(&String::from_utf8_lossy(v))).collect::<Vec<_>>().join(", "),
        trimmed.load(Relaxed),
        read_stats.records,
        read_stats.bases,
        wstats.kmers,
        wstats.tracked,
        bloom.bits(),
        inputs.iter().map(|p| json_str(p)).collect::<Vec<_>>().join(", ")
    );
    let db_path = format!("{out}.cpkdb");
    let mut w = BufWriter::with_capacity(8 << 20, File::create(&db_path).map_err(|e| format!("cannot write {db_path}: {e}"))?);
    w.write_all(DB_MAGIC)?;
    w.write_all(&(header.len() as u32).to_le_bytes())?;
    w.write_all(header.as_bytes())?;
    let group = threads * 2;
    let mut ci = 0;
    while ci < nchunks {
        let hi = (ci + group).min(nchunks);
        let bufs: Vec<Vec<u8>> = (ci..hi)
            .into_par_iter()
            .map(|c| {
                let mut b = Vec::with_capacity(CH * 4);
                for i in c * CH..((c + 1) * CH).min(slots) {
                    let key = table.keys[i].load(Relaxed);
                    if key != EMPTY {
                        let cnt = counts[i].load(Relaxed);
                        if cnt >= min_count {
                            b.extend_from_slice(&key.to_le_bytes());
                            b.extend_from_slice(&cnt.to_le_bytes());
                        }
                    }
                }
                b
            })
            .collect();
        for b in bufs {
            w.write_all(&b)?;
        }
        ci = hi;
    }
    w.flush()?;

    let mut hw = BufWriter::new(File::create(format!("{out}.hist.tsv"))?);
    writeln!(hw, "count\tn_kmers")?;
    for (c, &n) in hist.iter().enumerate() {
        if n > 0 {
            writeln!(hw, "{c}\t{n}")?;
        }
    }
    hw.flush()?;
    let mut jw = File::create(format!("{out}.json"))?;
    writeln!(
        jw,
        "{{\"db\": {}, \"header\": {header}, \"seconds\": {:.1}, \"peak_rss_gb\": {:.2}}}",
        json_str(&db_path),
        crate::util::elapsed(),
        peak_rss_bytes() as f64 / 1e9
    )?;
    log!("count: wrote {db_path} ({n_entries} entries); peak RSS {:.2} GB", peak_rss_bytes() as f64 / 1e9);
    Ok(())
}

/// Parsed database header.
pub struct DbHeader {
    pub k: usize,
    pub sample: u64,
    pub n_entries: u64,
    pub json: String,
}

/// Value of an unsigned integer field in a flat JSON object written by this program.
pub fn json_u64(json: &str, key: &str) -> Option<u64> {
    let pat = format!("\"{key}\":");
    let i = json.find(&pat)? + pat.len();
    let rest = json[i..].trim_start();
    let end = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
    rest[..end].parse().ok()
}

pub fn open_db(path: &str) -> Res<(DbHeader, BufReader<File>)> {
    let f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let flen = f.metadata()?.len();
    let mut r = BufReader::with_capacity(8 << 20, f);
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).map_err(|_| format!("{path}: not a cpk database"))?;
    if &magic != DB_MAGIC {
        return Err(format!("{path}: not a cpk database (bad magic)").into());
    }
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let hl = u32::from_le_bytes(lb) as usize;
    let mut hb = vec![0u8; hl];
    r.read_exact(&mut hb)?;
    let json = String::from_utf8(hb).map_err(|_| format!("{path}: corrupt header"))?;
    let need = |key: &str| json_u64(&json, key).ok_or_else(|| format!("{path}: header lacks '{key}'"));
    let h = DbHeader { k: need("k")? as usize, sample: need("sample")?, n_entries: need("n_entries")?, json };
    let expect = 12 + hl as u64 + h.n_entries * DB_ENTRY as u64;
    if flen != expect {
        return Err(format!("{path}: file is {flen} bytes but the header implies {expect}; the database is truncated or corrupt").into());
    }
    Ok((h, r))
}

/// `cpk dump --db X.cpkdb`: print "kmer<TAB>count" for every entry (testing and inspection).
pub fn dump(a: &Args) -> Res<()> {
    let path = a.req("db")?;
    let min = a.usize_or("min-count", 0)? as u16;
    let (h, mut r) = open_db(path)?;
    let stdout = std::io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut e = [0u8; DB_ENTRY];
    for _ in 0..h.n_entries {
        r.read_exact(&mut e)?;
        let hash = u64::from_le_bytes(e[0..8].try_into().unwrap());
        let c = u16::from_le_bytes(e[8..10].try_into().unwrap());
        if c >= min {
            writeln!(w, "{}\t{}", decode(fmix64_inv(hash), h.k), c)?;
        }
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finders() -> Vec<(memchr::memmem::Finder<'static>, Vec<u8>)> {
        vec![(memchr::memmem::Finder::new(TRUSEQ).into_owned(), TRUSEQ.to_vec())]
    }

    #[test]
    fn adapter_is_cut() {
        let f = finders();
        // a full adapter start inside the read
        assert_eq!(trim_adapter(b"ACGTACGTTTAGATCGGAAGAGCACACGTC", &f), b"ACGTACGTTT");
        // a partial adapter start (at least 5 bases) at the read's end
        assert_eq!(trim_adapter(b"ACGTACGTTTAGATCGG", &f), b"ACGTACGTTT");
        assert_eq!(trim_adapter(b"ACGTACGTTTAGATC", &f), b"ACGTACGTTT");
        // four bases are too few to call
        assert_eq!(trim_adapter(b"ACGTACGTTTAGAT", &f), b"ACGTACGTTTAGAT");
        // nothing to cut
        assert_eq!(trim_adapter(b"ACGTACGTTTACGT", &f), b"ACGTACGTTTACGT");
        // a read that is all adapter
        assert_eq!(trim_adapter(b"AGATCGGAAGAGCACACG", &f), b"");
        assert_eq!(trim_adapter(b"", &f), b"");
        assert_eq!(trim_adapter(b"ACGT", &[]), b"ACGT");
    }

    #[test]
    fn header_numbers() {
        let j = "{\"k\": 31, \"sample\": 4, \"n_entries\": 12345}";
        assert_eq!(json_u64(j, "k"), Some(31));
        assert_eq!(json_u64(j, "n_entries"), Some(12345));
        assert_eq!(json_u64(j, "absent"), None);
    }
}
