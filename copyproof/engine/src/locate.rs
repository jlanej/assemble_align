//! `cpk locate`: where do the k-mers of some query sequences occur in the release?
//!
//! Writes one line per (query k-mer, assembly occurrence): query, query position, sequence,
//! position, relative strand and the k-mer's total number of occurrences in the release.
//! k-mers with more than `--max-occ` occurrences are not listed: they cannot anchor a read.
//!
//! The queries are the sequences of a FASTA (`--query`) or stretches of the release itself
//! (`--regions`, a BED file with a name in column 4): the second form lists the other copies of
//! each stretch, and the stretch itself.

use crate::count::check_k_sample;
use crate::fastx::{load_named, Contig};
use crate::kmer::{fmix64, for_each_canon, EMPTY};
use crate::scan::{chunks_of, load_assemblies, tracked_kmers};
use crate::table::{atomic_u32_vec, KeyTable};
use crate::util::{labelled, threads_from, Args, Res};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::Ordering::Relaxed;

pub fn run(a: &Args) -> Res<()> {
    let k = a.usize_or("k", 31)?;
    let s = a.usize_or("sample", 1)? as u64;
    check_k_sample(k, s)?;
    let smask = s - 1;
    let threads = threads_from(a)?;
    let max_occ = a.usize_or("max-occ", 8)?.max(1) as u32;
    let out = a.req("out")?.to_string();
    let asm_specs = labelled(&a.all("asm"), "asm")?;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

    let asm = load_assemblies(&asm_specs, a.get("group-map"))?;
    let contigs = &asm.contigs;
    let queries = match (a.get("query"), a.get("regions")) {
        (Some(q), None) => load_named(q)?,
        (None, Some(r)) => regions_of(contigs, r)?,
        _ => return Err("locate: give --query Q.fa or --regions R.bed".into()),
    };
    let total: usize = queries.iter().map(|q| q.1.len().saturating_sub(k - 1)).sum();
    if total * 2 + 1024 > u32::MAX as usize {
        return Err("locate: too many query k-mers; pass fewer queries or a larger --sample".into());
    }
    let table = KeyTable::new(total * 2 / s as usize + 1024);
    queries.par_iter().for_each(|(_, seq)| {
        for_each_canon(seq, k, |_, canon, _| {
            let h = fmix64(canon);
            if h & smask == 0 && h != EMPTY {
                table.insert(h);
            }
        });
    });
    let occ = atomic_u32_vec(table.len());
    log!("locate: {} queries, {} k-mer positions", queries.len(), total);

    let chunks = chunks_of(contigs);
    // (slot, contig, position, forward strand is canonical)
    let mut hits: Vec<(u32, u32, u32, bool)> = chunks
        .par_iter()
        .map_init(Vec::<u8>::new, |codes, &(ci, st, en)| {
            let mut v = vec![];
            tracked_kmers(&contigs[ci], st, en, k, smask, codes, |p, h, fw| {
                if let Some(slot) = table.find(h) {
                    if occ[slot].fetch_add(1, Relaxed) < max_occ {
                        v.push((slot as u32, ci as u32, p as u32, fw));
                    }
                }
            });
            v
        })
        .flatten()
        .collect();
    hits.par_sort_unstable();
    log!("locate: {} listed occurrences in the release", hits.len());

    let mut w = BufWriter::with_capacity(4 << 20, File::create(&out).map_err(|e| format!("cannot write {out}: {e}"))?);
    writeln!(w, "query\tqpos\tcontig\tpos\tstrand\tocc")?;
    for grp in queries.chunks(4096) {
        let texts: Vec<String> = grp
            .par_iter()
            .map(|(name, seq)| {
                let mut t = String::new();
                for_each_canon(seq, k, |p, canon, qfw| {
                    let h = fmix64(canon);
                    if h & smask != 0 || h == EMPTY {
                        return;
                    }
                    let slot = match table.find(h) {
                        Some(x) => x as u32,
                        None => return,
                    };
                    let n = occ[slot as usize].load(Relaxed);
                    if n == 0 || n > max_occ {
                        return; // absent from the release, or too repeated to anchor anything
                    }
                    let lo = hits.partition_point(|x| x.0 < slot);
                    for x in hits[lo..].iter().take_while(|x| x.0 == slot) {
                        let strand = if x.3 == qfw { '+' } else { '-' };
                        t.push_str(&format!("{name}\t{p}\t{}\t{}\t{strand}\t{n}\n", contigs[x.1 as usize].name, x.2));
                    }
                });
                t
            })
            .collect();
        for t in texts {
            w.write_all(t.as_bytes())?;
        }
    }
    w.flush()?;
    log!("locate: wrote {out}");
    Ok(())
}

/// The stretches of the release named in a BED file (sequence, start, end, name), as ASCII.
fn regions_of(contigs: &[Contig], path: &str) -> Res<Vec<(String, Vec<u8>)>> {
    let index: HashMap<&str, usize> = contigs.iter().enumerate().map(|(i, c)| (c.name.as_str(), i)).collect();
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut out = vec![];
    let mut codes = vec![];
    for (ln, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        let bad = || format!("{path} line {}: expected sequence, start, end, name", ln + 1);
        if t.len() < 4 {
            return Err(bad().into());
        }
        let ci = *index.get(t[0]).ok_or_else(|| format!("{path} line {}: no sequence {} in the release", ln + 1, t[0]))?;
        let st: usize = t[1].parse().map_err(|_| bad())?;
        let en: usize = t[2].parse().map_err(|_| bad())?;
        contigs[ci].codes(st, en, &mut codes);
        out.push((t[3].to_string(), codes.iter().map(|&c| b"ACGTN"[c as usize]).collect()));
    }
    Ok(out)
}

/// `cpk extract`: the sequence of stretches of the release (BED with a name in column 4), as FASTA.
pub fn extract(a: &Args) -> Res<()> {
    let out = a.req("out")?.to_string();
    let asm_specs = labelled(&a.all("asm"), "asm")?;
    let asm = load_assemblies(&asm_specs, a.get("group-map"))?;
    let regions = regions_of(&asm.contigs, a.req("regions")?)?;
    let mut w = BufWriter::with_capacity(4 << 20, File::create(&out).map_err(|e| format!("cannot write {out}: {e}"))?);
    for (name, seq) in &regions {
        writeln!(w, ">{name}")?;
        w.write_all(seq)?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    log!("extract: wrote {} sequences to {out}", regions.len());
    Ok(())
}
