//! `cpk scan`: lay read k-mer counts along an assembly.
//!
//! For every tracked k-mer position of every sequence in the release, write the k-mer's
//! multiplicity in each assembly group (haplotype 1, haplotype 2, other) and its count in each read
//! set. Also list the read k-mers of the first read set that occur nowhere in the release.

use crate::count::{check_k_sample, open_db, DB_ENTRY};
use crate::fastx::{load_fasta, Contig};
use crate::kmer::{fmix64, fmix64_inv, for_each_canon_codes, EMPTY};
use crate::table::{atomic_u16_vec, atomic_u32_vec, Hll, KeyTable};
use crate::util::{json_str, labelled, peak_rss_bytes, threads_from, Args, Res};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::sync::atomic::{AtomicU16, AtomicU32, Ordering::Relaxed};

pub const CHUNK: usize = 1_000_000;
pub const BIN: usize = 100;
pub const MAX_GROUPS: usize = 4;

pub struct Assembly {
    pub contigs: Vec<Contig>,
    pub labels: Vec<String>,
}

/// Load every FASTA of a release. `specs` are (group label, path); `group_map` optionally assigns
/// individual sequences (by name) to a group, overriding their file's label.
pub fn load_assemblies(specs: &[(String, String)], group_map: Option<&str>) -> Res<Assembly> {
    if specs.is_empty() {
        return Err("give at least one --asm label=path".into());
    }
    let mut map: HashMap<String, String> = HashMap::new();
    if let Some(p) = group_map {
        let f = BufReader::new(File::open(p).map_err(|e| format!("cannot open {p}: {e}"))?);
        for line in f.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut it = line.split('\t');
            match (it.next(), it.next()) {
                (Some(n), Some(l)) => {
                    map.insert(n.to_string(), l.to_string());
                }
                _ => return Err(format!("{p}: expected 'sequence<TAB>group', got '{line}'").into()),
            }
        }
    }
    let loaded: Vec<Res<Vec<Contig>>> = specs.par_iter().map(|(_, path)| load_fasta(path, |_| Ok(0))).collect();
    let mut contigs: Vec<Contig> = vec![];
    let mut file_of: Vec<usize> = vec![];
    for (fi, r) in loaded.into_iter().enumerate() {
        let v = r?;
        file_of.extend(std::iter::repeat(fi).take(v.len()));
        contigs.extend(v);
    }
    // Groups are numbered in order of first appearance among the sequences.
    let mut labels: Vec<String> = vec![];
    for (c, &fi) in contigs.iter_mut().zip(file_of.iter()) {
        let label = map.get(&c.name).unwrap_or(&specs[fi].0);
        let gi = match labels.iter().position(|x| x == label) {
            Some(i) => i,
            None => {
                labels.push(label.clone());
                labels.len() - 1
            }
        };
        c.group = gi as u8;
    }
    if labels.len() > MAX_GROUPS {
        return Err(format!("at most {MAX_GROUPS} assembly groups are supported, got {}: {:?}", labels.len(), labels).into());
    }
    let mut seen: HashSet<&str> = HashSet::new();
    let dup = contigs.iter().any(|c| !seen.insert(c.name.as_str()));
    if dup {
        log!("sequence names repeat across files; prefixing each name with its file's label");
        for (c, &fi) in contigs.iter_mut().zip(file_of.iter()) {
            c.name = format!("{}#{}", specs[fi].0, c.name);
        }
        let mut seen: HashSet<&str> = HashSet::new();
        if let Some(c) = contigs.iter().find(|c| !seen.insert(c.name.as_str())) {
            return Err(format!("sequence name '{}' occurs twice within one group", c.name).into());
        }
    }
    for c in &contigs {
        if c.len > u32::MAX as usize {
            return Err(format!("sequence {} is longer than 4.29 Gb, which is not supported", c.name).into());
        }
    }
    Ok(Assembly { contigs, labels })
}

/// (contig index, start, end) pieces covering every base, CHUNK bases each.
pub fn chunks_of(contigs: &[Contig]) -> Vec<(usize, usize, usize)> {
    let mut v = vec![];
    for (ci, c) in contigs.iter().enumerate() {
        let mut st = 0;
        while st < c.len {
            let en = (st + CHUNK).min(c.len);
            v.push((ci, st, en));
            st = en;
        }
    }
    v
}

/// Call `f(start, hash, forward_is_canonical)` for each tracked k-mer starting in [st, en).
#[inline]
pub fn tracked_kmers<F: FnMut(usize, u64, bool)>(c: &Contig, st: usize, en: usize, k: usize, smask: u64, codes: &mut Vec<u8>, mut f: F) {
    c.codes(st, (en + k - 1).min(c.len), codes);
    for_each_canon_codes(codes, k, |p, canon, fw| {
        let h = fmix64(canon);
        if h & smask == 0 && h != EMPTY {
            f(st + p, h, fw);
        }
    });
}

#[inline]
fn bump_group(m: &AtomicU32, g: u8) {
    let sh = 8 * g as u32;
    let mut cur = m.load(Relaxed);
    loop {
        if (cur >> sh) & 0xff == 0xff {
            return;
        }
        match m.compare_exchange_weak(cur, cur + (1 << sh), Relaxed, Relaxed) {
            Ok(_) => return,
            Err(x) => cur = x,
        }
    }
}

struct ReadSet {
    label: String,
    path: String,
    header_json: String,
    entries: u64,
    found: u64,
    absent: u64,
}

pub fn run(a: &Args) -> Res<()> {
    let k = a.usize_or("k", 31)?;
    let s = a.usize_or("sample", 4)? as u64;
    check_k_sample(k, s)?;
    let smask = s - 1;
    let threads = threads_from(a)?;
    let asm_specs = labelled(&a.all("asm"), "asm")?;
    let read_specs = labelled(&a.all("reads"), "reads")?;
    let out = a.req("out")?.to_string();
    let novel_min = a.usize_or("novel-min", 3)?.clamp(1, 65_535) as u16;
    let load = a.f64_or("load", 0.6)?.clamp(0.2, 0.85);
    let high_copy = a.usize_or("high-copy", 8)? as u32;
    let high_depth = a.f64_or("high-depth", 20.0)?;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

    // Check the read databases before the expensive work.
    for (label, path) in &read_specs {
        let (h, _) = open_db(path)?;
        if h.k != k {
            return Err(format!("read set '{label}' was counted with k={} but scan uses k={k}", h.k).into());
        }
        if h.sample > s {
            return Err(format!("read set '{label}' tracks 1/{} of k-mers, coarser than the requested 1/{s}", h.sample).into());
        }
    }

    log!("scan: loading {} FASTA file(s)", asm_specs.len());
    let asm = load_assemblies(&asm_specs, a.get("group-map"))?;
    let contigs = &asm.contigs;
    let ng = asm.labels.len();
    let nr = read_specs.len();
    let total_bases: usize = contigs.iter().map(|c| c.len).sum();
    log!("scan: {} sequences, {:.3} Gb, groups {:?}", contigs.len(), total_bases as f64 / 1e9, asm.labels);
    let chunks = chunks_of(contigs);

    // Distinct tracked k-mers, estimated, to size the table.
    let (hll, instances) = chunks
        .par_iter()
        .fold(
            || (Hll::new(), 0u64, Vec::<u8>::new()),
            |(mut h, mut n, mut codes), &(ci, st, en)| {
                tracked_kmers(&contigs[ci], st, en, k, smask, &mut codes, |_, hash, _| {
                    h.add(hash);
                    n += 1;
                });
                (h, n, codes)
            },
        )
        .map(|(h, n, _)| (h, n))
        .reduce(
            || (Hll::new(), 0u64),
            |mut x, y| {
                x.0.merge(&y.0);
                (x.0, x.1 + y.1)
            },
        );
    let est = hll.estimate();
    let slots = ((est * 1.03 / load) as usize + 4096).min(instances as usize * 2 + 4096);
    log!(
        "scan: {} tracked k-mer positions, about {:.0} distinct; table {:.2} GB",
        instances,
        est,
        slots as f64 * (12.0 + 2.0 * nr as f64) / 1e9
    );
    let table = KeyTable::new(slots);
    let mult = atomic_u32_vec(slots);
    let rcounts: Vec<Vec<AtomicU16>> = (0..nr).map(|_| atomic_u16_vec(slots)).collect();

    // Insert assembly k-mers with their multiplicity per group.
    let distinct: u64 = chunks
        .par_iter()
        .map_init(Vec::<u8>::new, |codes, &(ci, st, en)| {
            let c = &contigs[ci];
            let mut newly = 0u64;
            tracked_kmers(c, st, en, k, smask, codes, |_, hash, _| {
                let (slot, new) = table.insert(hash);
                newly += new as u64;
                bump_group(&mult[slot], c.group);
            });
            newly
        })
        .sum();
    log!("scan: {distinct} distinct tracked k-mers in the release (table load {:.2})", distinct as f64 / slots as f64);

    // Join the read databases.
    let mut sets: Vec<ReadSet> = vec![];
    let mut novel: Vec<(u64, u16)> = vec![];
    let mut novel_table: Option<KeyTable> = None;
    let mut novel_counts: Vec<Vec<AtomicU16>> = vec![];
    const BLOCK: usize = 4 << 20; // entries per read
    for (r, (label, path)) in read_specs.iter().enumerate() {
        let (h, mut rdr) = open_db(path)?;
        let (mut found, mut absent) = (0u64, 0u64);
        let mut left = h.n_entries as usize;
        let mut buf = vec![0u8; BLOCK * DB_ENTRY];
        while left > 0 {
            let n = left.min(BLOCK);
            rdr.read_exact(&mut buf[..n * DB_ENTRY])?;
            left -= n;
            let parts: Vec<(u64, u64, Vec<(u64, u16)>)> = buf[..n * DB_ENTRY]
                .par_chunks(DB_ENTRY * 32_768)
                .map(|ch| {
                    let (mut f, mut ab) = (0u64, 0u64);
                    let mut nov = vec![];
                    for e in ch.chunks_exact(DB_ENTRY) {
                        let hash = u64::from_le_bytes(e[0..8].try_into().unwrap());
                        let c = u16::from_le_bytes(e[8..10].try_into().unwrap());
                        if hash & smask != 0 {
                            continue;
                        }
                        match table.find(hash) {
                            Some(i) => {
                                rcounts[r][i].store(c, Relaxed);
                                f += 1;
                            }
                            None => {
                                ab += 1;
                                if r == 0 {
                                    if c >= novel_min {
                                        nov.push((hash, c));
                                    }
                                } else if let Some(nt) = &novel_table {
                                    if let Some(j) = nt.find(hash) {
                                        novel_counts[r - 1][j].store(c, Relaxed);
                                    }
                                }
                            }
                        }
                    }
                    (f, ab, nov)
                })
                .collect();
            for (f, ab, nov) in parts {
                found += f;
                absent += ab;
                novel.extend(nov);
            }
        }
        log!("scan: read set '{label}': {} entries, {found} in the release, {absent} absent", h.n_entries);
        sets.push(ReadSet { label: label.clone(), path: path.clone(), header_json: h.json.clone(), entries: h.n_entries, found, absent });
        if r == 0 {
            novel.par_sort_unstable();
            let nt = KeyTable::new(novel.len() * 2 + 1024);
            novel.par_iter().for_each(|&(hash, _)| {
                nt.insert(hash);
            });
            novel_counts = (1..nr).map(|_| atomic_u16_vec(nt.len())).collect();
            novel_table = Some(nt);
            log!("scan: {} tracked k-mers with count >= {novel_min} in '{label}' are absent from the release", novel.len());
        }
    }

    // Spectrum: distinct k-mers by multiplicity pattern (each group capped at 3) and first-set count.
    const CCAP: usize = 1024;
    let npat = 1usize << (2 * ng);
    let spectrum: Vec<u64> = (0..slots.div_ceil(1 << 22))
        .into_par_iter()
        .map(|ch| {
            let mut hst = vec![0u64; npat * CCAP];
            for i in ch << 22..((ch + 1) << 22).min(slots) {
                if table.keys[i].load(Relaxed) == EMPTY {
                    continue;
                }
                let m = mult[i].load(Relaxed);
                let mut pat = 0usize;
                for g in 0..ng {
                    pat |= (((m >> (8 * g)) & 0xff).min(3) as usize) << (2 * g);
                }
                let c = if nr > 0 { (rcounts[0][i].load(Relaxed) as usize).min(CCAP - 1) } else { 0 };
                hst[pat * CCAP + c] += 1;
            }
            hst
        })
        .reduce(
            || vec![0u64; npat * CCAP],
            |mut x, y| {
                for (p, q) in x.iter_mut().zip(y.iter()) {
                    *p += *q;
                }
                x
            },
        );
    {
        let mut w = BufWriter::new(File::create(format!("{out}.spectrum.tsv"))?);
        writeln!(w, "{}\tcount\tn_kmers", asm.labels.iter().map(|l| format!("m_{l}")).collect::<Vec<_>>().join("\t"))?;
        for pat in 0..npat {
            for c in 0..CCAP {
                let n = spectrum[pat * CCAP + c];
                if n > 0 {
                    let ms: Vec<String> = (0..ng).map(|g| ((pat >> (2 * g)) & 3).to_string()).collect();
                    writeln!(w, "{}\t{c}\t{n}", ms.join("\t"))?;
                }
            }
        }
    }

    // Hashes of tracked k-mers that are high-copy: above the multiplicity cap in the release, or
    // read far more often than the release's copies explain (mitochondrion, rDNA, collapsed
    // repeats). Around such a k-mer every single-base error is seen often enough to look like
    // solid sequence that the release lacks; `recruit` uses this set to ignore those.
    {
        // depth per copy: the median read count of the commonest multiplicity pattern
        let mut per_pat = vec![0u64; npat];
        for pat in 0..npat {
            per_pat[pat] = spectrum[pat * CCAP..(pat + 1) * CCAP].iter().sum();
        }
        let top = (0..npat).max_by_key(|&p| per_pat[p]).unwrap_or(0);
        let copies: u32 = (0..ng).map(|g| ((top >> (2 * g)) & 3) as u32).sum::<u32>().max(1);
        let nonzero: u64 = spectrum[top * CCAP + 1..(top + 1) * CCAP].iter().sum();
        let (mut acc, mut median) = (0u64, 0usize);
        for c in 1..CCAP {
            acc += spectrum[top * CCAP + c];
            if acc * 2 >= nonzero {
                median = c;
                break;
            }
        }
        let lambda = median as f64 / copies as f64;
        let high_count = if nr > 0 && lambda > 0.0 { (high_depth * lambda).ceil().min(65_535.0) as u16 } else { u16::MAX };
        let hashes: Vec<u64> = (0..slots.div_ceil(1 << 22))
            .into_par_iter()
            .flat_map_iter(|ch| {
                let mut v = vec![];
                for i in ch << 22..((ch + 1) << 22).min(slots) {
                    let key = table.keys[i].load(Relaxed);
                    if key == EMPTY {
                        continue;
                    }
                    let m = mult[i].load(Relaxed);
                    let tot: u32 = (0..ng).map(|g| (m >> (8 * g)) & 0xff).sum();
                    let c = if nr > 0 { rcounts[0][i].load(Relaxed) } else { 0 };
                    if tot > high_copy || c >= high_count {
                        v.push(key);
                    }
                }
                v
            })
            .collect();
        let mut w = BufWriter::with_capacity(4 << 20, File::create(format!("{out}.repeat.u64"))?);
        for h in &hashes {
            w.write_all(&h.to_le_bytes())?;
        }
        w.flush()?;
        log!(
            "scan: {} tracked k-mers are high-copy (more than {high_copy} in the release, or read count >= {high_count}; depth per copy {lambda:.1})",
            hashes.len()
        );
    }

    // Positional track.
    let rec = 1 + ng + 2 * nr;
    let pay = rec - 1;
    let track_path = format!("{out}.track.bin");
    let mut tw = BufWriter::with_capacity(16 << 20, File::create(&track_path).map_err(|e| format!("cannot write {track_path}: {e}"))?);
    let theader = format!(
        "{{\"format\": \"cpktrack\", \"version\": 1, \"k\": {k}, \"sample\": {s}, \"groups\": [{}], \"reads\": [{}], \"record_bytes\": {rec}, \"bin\": {BIN}}}",
        asm.labels.iter().map(|l| json_str(l)).collect::<Vec<_>>().join(", "),
        read_specs.iter().map(|(l, _)| json_str(l)).collect::<Vec<_>>().join(", ")
    );
    tw.write_all(b"CPKTRK01")?;
    tw.write_all(&(theader.len() as u32).to_le_bytes())?;
    tw.write_all(theader.as_bytes())?;
    let mut gcw = BufWriter::with_capacity(4 << 20, File::create(format!("{out}.gc{BIN}.u8"))?);
    let mut acw = BufWriter::with_capacity(4 << 20, File::create(format!("{out}.acgt{BIN}.u8"))?);
    let mut cw = BufWriter::new(File::create(format!("{out}.contigs.tsv"))?);
    writeln!(cw, "idx\tname\tgroup\tlength\trec_start\tn_rec\tn_kmers\tbin_start\tn_bins")?;
    let mut gw = BufWriter::new(File::create(format!("{out}.gaps.bed"))?);

    let mut rec_index = 0u64;
    let mut bin_index = 0u64;
    let filler = {
        let mut f = vec![0u8; rec];
        f[0] = 255;
        f
    };
    for (ci, c) in contigs.iter().enumerate() {
        for &(s0, e0) in &c.nruns {
            writeln!(gw, "{}\t{s0}\t{e0}", c.name)?;
        }
        let rec_start = rec_index;
        let bin_start = bin_index;
        let mut n_kmers = 0u64;
        let mut prev = 0usize;
        let cchunks: Vec<(usize, usize)> = (0..c.len.div_ceil(CHUNK)).map(|i| (i * CHUNK, ((i + 1) * CHUNK).min(c.len))).collect();
        for grp in cchunks.chunks(threads * 4) {
            let parts: Vec<(Vec<u32>, Vec<u8>, Vec<u8>, Vec<u8>)> = grp
                .par_iter()
                .map_init(Vec::<u8>::new, |codes, &(st, en)| {
                    let mut pos: Vec<u32> = Vec::with_capacity((en - st) / s as usize + 64);
                    let mut payload: Vec<u8> = Vec::with_capacity(((en - st) / s as usize + 64) * pay);
                    tracked_kmers(c, st, en, k, smask, codes, |p, hash, _| {
                        let slot = table.find(hash).expect("assembly k-mer missing from its own table");
                        let m = mult[slot].load(Relaxed);
                        pos.push((p - st) as u32);
                        for g in 0..ng {
                            payload.push(((m >> (8 * g)) & 0xff) as u8);
                        }
                        for rc in rcounts.iter() {
                            payload.extend_from_slice(&rc[slot].load(Relaxed).to_le_bytes());
                        }
                    });
                    // base composition per bin over [st, en); `codes` starts at st
                    let nb = (en - st).div_ceil(BIN);
                    let mut gc = vec![0u8; nb];
                    let mut acgt = vec![0u8; nb];
                    for (i, &code) in codes[..en - st].iter().enumerate() {
                        if code < 4 {
                            acgt[i / BIN] += 1;
                            gc[i / BIN] += (code == 1 || code == 2) as u8;
                        }
                    }
                    (pos, payload, gc, acgt)
                })
                .collect();
            for (i, (pos, payload, gc, acgt)) in parts.iter().enumerate() {
                let base = grp[i].0;
                for (j, &p) in pos.iter().enumerate() {
                    let abs = base + p as usize;
                    let mut gap = abs - prev;
                    while gap > 255 {
                        tw.write_all(&filler)?;
                        rec_index += 1;
                        gap -= 255;
                    }
                    tw.write_all(&[gap as u8])?;
                    tw.write_all(&payload[j * pay..(j + 1) * pay])?;
                    rec_index += 1;
                    prev = abs;
                }
                n_kmers += pos.len() as u64;
                gcw.write_all(gc)?;
                acw.write_all(acgt)?;
                bin_index += gc.len() as u64;
            }
        }
        writeln!(
            cw,
            "{ci}\t{}\t{}\t{}\t{rec_start}\t{}\t{n_kmers}\t{bin_start}\t{}",
            c.name,
            asm.labels[c.group as usize],
            c.len,
            rec_index - rec_start,
            bin_index - bin_start
        )?;
    }
    tw.flush()?;
    gcw.flush()?;
    acw.flush()?;
    cw.flush()?;
    gw.flush()?;
    log!("scan: wrote {track_path} ({rec_index} records of {rec} bytes)");

    // Read k-mers absent from the release.
    if nr > 0 {
        let npath = format!("{out}.novel.bin");
        let mut nw = BufWriter::with_capacity(8 << 20, File::create(&npath)?);
        let nheader = format!(
            "{{\"format\": \"cpknovel\", \"version\": 1, \"k\": {k}, \"sample\": {s}, \"n\": {}, \"min_count\": {novel_min}, \"reads\": [{}], \"record_bytes\": {}}}",
            novel.len(),
            read_specs.iter().map(|(l, _)| json_str(l)).collect::<Vec<_>>().join(", "),
            8 + 2 * nr
        );
        nw.write_all(b"CPKNOV01")?;
        nw.write_all(&(nheader.len() as u32).to_le_bytes())?;
        nw.write_all(nheader.as_bytes())?;
        let nt = novel_table.as_ref().unwrap();
        for &(hash, c) in &novel {
            nw.write_all(&fmix64_inv(hash).to_le_bytes())?;
            nw.write_all(&c.to_le_bytes())?;
            if nr > 1 {
                let j = nt.find(hash).unwrap();
                for pc in novel_counts.iter() {
                    nw.write_all(&pc[j].load(Relaxed).to_le_bytes())?;
                }
            }
        }
        nw.flush()?;
    }

    // Summary.
    let mut groups_json = vec![];
    for (g, l) in asm.labels.iter().enumerate() {
        let (mut n, mut b, mut nn) = (0u64, 0u64, 0u64);
        for c in contigs.iter().filter(|c| c.group as usize == g) {
            n += 1;
            b += c.len as u64;
            nn += c.nruns.iter().map(|r| (r.1 - r.0) as u64).sum::<u64>();
        }
        groups_json.push(format!("{{\"label\": {}, \"sequences\": {n}, \"bases\": {b}, \"non_acgt\": {nn}}}", json_str(l)));
    }
    let sets_json: Vec<String> = sets
        .iter()
        .map(|x| {
            format!(
                "{{\"label\": {}, \"path\": {}, \"entries\": {}, \"in_release\": {}, \"absent\": {}, \"db\": {}}}",
                json_str(&x.label),
                json_str(&x.path),
                x.entries,
                x.found,
                x.absent,
                x.header_json
            )
        })
        .collect();
    let mut jw = File::create(format!("{out}.scan.json"))?;
    writeln!(
        jw,
        "{{\"k\": {k}, \"sample\": {s}, \"asm\": [{}], \"groups\": [{}], \"reads\": [{}], \"tracked_positions\": {instances}, \"distinct_tracked\": {distinct}, \"hll_estimate\": {:.0}, \"table_slots\": {slots}, \"track_records\": {rec_index}, \"record_bytes\": {rec}, \"novel_min\": {novel_min}, \"n_novel\": {}, \"seconds\": {:.1}, \"peak_rss_gb\": {:.2}}}",
        asm_specs.iter().map(|(l, p)| format!("{{\"label\": {}, \"path\": {}}}", json_str(l), json_str(p))).collect::<Vec<_>>().join(", "),
        groups_json.join(", "),
        sets_json.join(", "),
        est,
        novel.len(),
        crate::util::elapsed(),
        peak_rss_bytes() as f64 / 1e9
    )?;
    log!("scan: done; peak RSS {:.2} GB", peak_rss_bytes() as f64 / 1e9);
    Ok(())
}
