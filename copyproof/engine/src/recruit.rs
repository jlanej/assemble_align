//! `cpk recruit`: one streaming pass over a read set, keeping reads that carry k-mers of a given set
//! (for example solid read k-mers that are absent from the assembly).

use crate::count::check_k_sample;
use crate::fastx::{stream_reads, Batch};
use crate::kmer::{fmix64, for_each_canon, EMPTY};
use crate::table::KeyTable;
use crate::util::{threads_from, Args, Res};
use crossbeam_channel::bounded;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::sync::Mutex;

/// Read a file of little-endian u64 canonical k-mers.
pub fn load_kmer_file(path: &str) -> Res<Vec<u64>> {
    let mut f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let mut bytes = vec![];
    f.read_to_end(&mut bytes)?;
    if bytes.len() % 8 != 0 {
        return Err(format!("{path}: length is not a multiple of 8 bytes").into());
    }
    Ok(bytes.chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().unwrap())).collect())
}

pub fn run(a: &Args) -> Res<()> {
    let k = a.usize_or("k", 31)?;
    check_k_sample(k, 1)?;
    let threads = threads_from(a)?;
    let min_hits = a.usize_or("min-hits", 1)?.max(1);
    let max_reads = a.usize_or("max-reads", usize::MAX)? as u64;
    // Optional: hashes of tracked k-mers that are high-copy in the release (scan's .repeat.u64).
    // A read whose tracked k-mers are mostly of this kind is not kept.
    let sample = a.usize_or("sample", 4)? as u64;
    let max_repeat = a.f64_or("max-repeat-frac", 0.5)?;
    let repeat: Option<KeyTable> = match a.get("repeat-hashes") {
        Some(p) => {
            let hs = load_kmer_file(p)?;
            let t = KeyTable::new(hs.len() * 2 + 1024);
            for h in hs {
                if h != EMPTY {
                    t.insert(h);
                }
            }
            Some(t)
        }
        None => None,
    };
    let smask = sample.max(1) - 1;
    // If every k-mer of the set passes the tracking test of this rate, other k-mers need no lookup.
    let set_sample = a.usize_or("set-sample", 1)? as u64;
    if set_sample == 0 || !set_sample.is_power_of_two() {
        return Err("--set-sample must be a power of two".into());
    }
    let setmask = set_sample - 1;
    let out = a.req("out")?.to_string();
    let kmers = load_kmer_file(a.req("kmers")?)?;
    let inputs = a.pos.clone();
    if inputs.is_empty() {
        return Err("recruit: give at least one FASTQ/FASTA input ('-' for stdin)".into());
    }
    if kmers.len() > u32::MAX as usize {
        return Err("recruit: too many k-mers in the set".into());
    }
    if let Some(bad) = kmers.iter().find(|&&km| fmix64(km) & setmask != 0) {
        return Err(format!("--set-sample {set_sample} does not hold for k-mer {bad} of the set").into());
    }
    let table = KeyTable::new(kmers.len() * 2 + 1024);
    let mut index = vec![0u32; table.len()];
    for (i, &km) in kmers.iter().enumerate() {
        let h = fmix64(km);
        if h == EMPTY {
            continue;
        }
        let (slot, _) = table.insert(h);
        index[slot] = i as u32;
    }
    log!("recruit: {} k-mers in the set; k={k}; min hits {min_hits}", kmers.len());

    let writer = Mutex::new(BufWriter::with_capacity(4 << 20, File::create(&out).map_err(|e| format!("cannot write {out}: {e}"))?));
    let kept = AtomicU64::new(0);
    let seen = AtomicU64::new(0);
    let in_repeats = AtomicU64::new(0);
    let finders: Vec<(memchr::memmem::Finder<'static>, Vec<u8>)> = crate::count::adapter_finders(a, false);
    let nreaders = inputs.len().min((threads / 2).max(1));
    let next_input = AtomicUsize::new(0);
    let (tx, rx) = bounded::<Batch>(threads * 4);
    let mut reader_err: Option<String> = None;
    std::thread::scope(|sc| {
        let readers: Vec<_> = (0..nreaders)
            .map(|_| {
                let tx = tx.clone();
                let (inputs, next_input) = (&inputs, &next_input);
                sc.spawn(move || -> Result<(), String> {
                    loop {
                        let i = next_input.fetch_add(1, Relaxed);
                        if i >= inputs.len() {
                            return Ok(());
                        }
                        stream_reads(&inputs[i], i, 4 << 20, 0, &tx).map_err(|e| e.to_string())?;
                    }
                })
            })
            .collect();
        drop(tx);
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let rx = rx.clone();
                let (table, index, writer, kept, seen) = (&table, &index, &writer, &kept, &seen);
                let (repeat, in_repeats, finders) = (&repeat, &in_repeats, &finders);
                sc.spawn(move || {
                    let mut obuf: Vec<u8> = Vec::with_capacity(1 << 20);
                    let mut hits: Vec<(u32, u32)> = vec![];
                    let mut reps: Vec<u32> = vec![];
                    for batch in rx.iter() {
                        if kept.load(Relaxed) >= max_reads {
                            continue;
                        }
                        let mut nrec = 0u64;
                        let mut nkept = 0u64;
                        let mut nrep = 0u64;
                        for (ri, seq) in batch.data.split(|&b| b == b'\n').enumerate() {
                            if seq.is_empty() {
                                continue;
                            }
                            nrec += 1;
                            hits.clear();
                            let seq = crate::count::trim_adapter(seq, finders);
                            let (mut tracked, mut rep) = (0u32, 0u32);
                            reps.clear();
                            for_each_canon(seq, k, |p, canon, _| {
                                let h = fmix64(canon);
                                // with --set-sample the set holds tracked k-mers only: skip the rest unseen
                                if h & setmask == 0 {
                                    if let Some(slot) = table.find(h) {
                                        hits.push((index[slot], p as u32));
                                    }
                                }
                                if h & smask == 0 {
                                    tracked += 1;
                                    if let Some(r) = repeat {
                                        if r.find(h).is_some() {
                                            rep += 1;
                                            reps.push(p as u32);
                                        }
                                    }
                                }
                            });
                            let before = hits.len();
                            if before >= min_hits && rep as f64 > max_repeat * tracked as f64 {
                                nrep += 1;
                                continue;
                            }
                            // A set k-mer that sits among high-copy k-mers is most likely a sequencing
                            // error of one of them: every neighbour of a k-mer that occurs thousands of
                            // times is seen often enough to look solid.
                            if !reps.is_empty() {
                                let reach = 2 * k as u32;
                                hits.retain(|&(_, p)| !reps.iter().any(|&q| q.abs_diff(p) <= reach));
                            }
                            if before >= min_hits && hits.len() < min_hits {
                                nrep += 1;
                            }
                            if hits.len() >= min_hits {
                                nkept += 1;
                                write!(obuf, ">{}_{} hits={} k=", batch.file, batch.first_record + ri as u64, hits.len()).unwrap();
                                for (j, (ki, p)) in hits.iter().enumerate() {
                                    if j > 0 {
                                        obuf.push(b',');
                                    }
                                    write!(obuf, "{ki}:{p}").unwrap();
                                }
                                obuf.push(b'\n');
                                obuf.extend_from_slice(seq);
                                obuf.push(b'\n');
                            }
                        }
                        seen.fetch_add(nrec, Relaxed);
                        kept.fetch_add(nkept, Relaxed);
                        in_repeats.fetch_add(nrep, Relaxed);
                        if !obuf.is_empty() {
                            writer.lock().unwrap().write_all(&obuf).unwrap();
                            obuf.clear();
                        }
                    }
                })
            })
            .collect();
        drop(rx);
        for r in readers {
            if let Err(e) = r.join().unwrap() {
                reader_err = Some(e);
            }
        }
        for w in workers {
            w.join().unwrap();
        }
    });
    if let Some(e) = reader_err {
        return Err(e.into());
    }
    writer.lock().unwrap().flush()?;
    log!(
        "recruit: kept {} of {} reads; {} more carried set k-mers only inside high-copy repeat",
        kept.load(Relaxed),
        seen.load(Relaxed),
        in_repeats.load(Relaxed)
    );
    if kept.load(Relaxed) >= max_reads {
        log!("recruit: stopped keeping reads at --max-reads {max_reads}");
    }
    Ok(())
}
