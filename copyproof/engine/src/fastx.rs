//! Sequence input: whole assemblies held 2-bit packed in memory, and reads streamed as batches.

use crate::kmer::ENC;
use crate::util::Res;
use crossbeam_channel::Sender;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

/// Open a file (or stdin for "-"), transparently decompressing gzip/bgzip.
pub fn open_maybe_gz(path: &str) -> Res<Box<dyn BufRead + Send>> {
    let raw: Box<dyn Read + Send> = if path == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?)
    };
    let mut br = BufReader::with_capacity(1 << 20, raw);
    let head = br.fill_buf()?;
    if head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        Ok(Box::new(BufReader::with_capacity(4 << 20, MultiGzDecoder::new(br))))
    } else {
        Ok(Box::new(br))
    }
}

/// One assembly sequence, 2 bits per base, with runs of non-ACGT bases kept as intervals.
pub struct Contig {
    pub name: String,
    pub group: u8,
    pub len: usize,
    pub packed: Vec<u8>,
    /// Half-open runs of bases that are not A/C/G/T, sorted.
    pub nruns: Vec<(usize, usize)>,
}

impl Contig {
    /// Fill `out` with the codes (0..=3; 4 for non-ACGT) of bases [start, end).
    pub fn codes(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.clear();
        let end = end.min(self.len);
        if start >= end {
            return;
        }
        out.reserve(end - start);
        for i in start..end {
            out.push((self.packed[i >> 2] >> ((i & 3) * 2)) & 3);
        }
        // overlay non-ACGT runs
        let first = self.nruns.partition_point(|r| r.1 <= start);
        for &(s, e) in &self.nruns[first..] {
            if s >= end {
                break;
            }
            let a = s.max(start) - start;
            let b = e.min(end) - start;
            for x in &mut out[a..b] {
                *x = 4;
            }
        }
    }
}

/// Read every sequence of a FASTA file. `group_of` maps a sequence name to its group index.
pub fn load_fasta<F: Fn(&str) -> Res<u8>>(path: &str, group_of: F) -> Res<Vec<Contig>> {
    let mut rdr = open_maybe_gz(path)?;
    let mut contigs: Vec<Contig> = vec![];
    let mut line: Vec<u8> = Vec::with_capacity(1 << 16);

    struct Cur {
        name: String,
        group: u8,
        len: usize,
        packed: Vec<u8>,
        nruns: Vec<(usize, usize)>,
        byte: u8,
        nstart: Option<usize>,
    }
    fn finish(c: Cur, out: &mut Vec<Contig>) {
        let mut c = c;
        if c.len & 3 != 0 {
            c.packed.push(c.byte);
        }
        if let Some(s) = c.nstart {
            c.nruns.push((s, c.len));
        }
        c.packed.shrink_to_fit();
        out.push(Contig { name: c.name, group: c.group, len: c.len, packed: c.packed, nruns: c.nruns });
    }

    let mut cur: Option<Cur> = None;
    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            if let Some(c) = cur.take() {
                finish(c, &mut contigs);
            }
            let text = String::from_utf8_lossy(&line[1..]);
            let name = text.split_whitespace().next().unwrap_or("").to_string();
            if name.is_empty() {
                return Err(format!("{path}: FASTA record with an empty name").into());
            }
            let group = group_of(&name)?;
            cur = Some(Cur { name, group, len: 0, packed: vec![], nruns: vec![], byte: 0, nstart: None });
            continue;
        }
        let c = match cur.as_mut() {
            Some(c) => c,
            None => {
                if line.iter().all(|b| b.is_ascii_whitespace()) {
                    continue;
                }
                return Err(format!("{path}: sequence data before the first '>' header").into());
            }
        };
        for &b in line.iter() {
            let code = ENC[b as usize];
            if code == 4 {
                if b.is_ascii_whitespace() {
                    continue;
                }
                if c.nstart.is_none() {
                    c.nstart = Some(c.len);
                }
            } else if let Some(s) = c.nstart.take() {
                c.nruns.push((s, c.len));
            }
            c.byte |= (code & 3) << ((c.len & 3) * 2);
            c.len += 1;
            if c.len & 3 == 0 {
                c.packed.push(c.byte);
                c.byte = 0;
            }
        }
    }
    if let Some(c) = cur.take() {
        finish(c, &mut contigs);
    }
    if contigs.is_empty() {
        return Err(format!("{path}: no FASTA records").into());
    }
    Ok(contigs)
}

/// Sequences of consecutive records from one input, each terminated by '\n'.
pub struct Batch {
    pub file: usize,
    /// 0-based index, within the file, of the first record in this batch.
    pub first_record: u64,
    pub data: Vec<u8>,
}

#[derive(Default, Clone, Copy)]
pub struct ReadStats {
    pub records: u64,
    pub bases: u64,
}

/// Stream one FASTQ or FASTA input (optionally gzipped; "-" is stdin) as batches of sequences.
///
/// With `overlap` > 0, a FASTA record longer than a batch (a chromosome) is cut into pieces that
/// overlap by `overlap` bases; with overlap = k - 1 every k-mer falls in exactly one piece.
pub fn stream_reads(path: &str, file: usize, batch_bytes: usize, overlap: usize, tx: &Sender<Batch>) -> Res<ReadStats> {
    let mut rdr = open_maybe_gz(path)?;
    let first = {
        let head = rdr.fill_buf()?;
        if head.is_empty() {
            return Ok(ReadStats::default());
        }
        head[0]
    };
    let mut stats = ReadStats::default();
    let mut data: Vec<u8> = Vec::with_capacity(batch_bytes + (1 << 16));
    let mut first_record = 0u64;
    let mut line: Vec<u8> = Vec::with_capacity(1 << 12);

    macro_rules! flush {
        () => {
            if !data.is_empty() {
                let full = std::mem::replace(&mut data, Vec::with_capacity(batch_bytes + (1 << 16)));
                if tx.send(Batch { file, first_record, data: full }).is_err() {
                    return Err("worker threads stopped early".into());
                }
                #[allow(unused_assignments)]
                {
                    first_record = stats.records;
                }
            }
        };
    }

    match first {
        b'@' => {
            let mut li = 0u8;
            loop {
                line.clear();
                if rdr.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                match li {
                    0 => {
                        if line[0] != b'@' {
                            return Err(format!(
                                "{path}: record {} does not start with '@' (only 4-line FASTQ is supported)",
                                stats.records + 1
                            )
                            .into());
                        }
                    }
                    1 => {
                        let mut n = line.len();
                        while n > 0 && (line[n - 1] == b'\n' || line[n - 1] == b'\r') {
                            n -= 1;
                        }
                        data.extend_from_slice(&line[..n]);
                        data.push(b'\n');
                        stats.records += 1;
                        stats.bases += n as u64;
                    }
                    2 => {
                        if line[0] != b'+' {
                            return Err(format!(
                                "{path}: record {} has no '+' line (only 4-line FASTQ is supported)",
                                stats.records
                            )
                            .into());
                        }
                    }
                    _ => {
                        if data.len() >= batch_bytes {
                            flush!();
                        }
                    }
                }
                li = (li + 1) & 3;
            }
            if li != 0 {
                return Err(format!("{path}: truncated FASTQ (the last record is incomplete)").into());
            }
        }
        b'>' => {
            let mut open = false;
            loop {
                line.clear();
                if rdr.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                if line[0] == b'>' {
                    if open {
                        data.push(b'\n');
                        if data.len() >= batch_bytes {
                            flush!();
                        }
                    }
                    open = true;
                    stats.records += 1;
                    continue;
                }
                let mut n = line.len();
                while n > 0 && (line[n - 1] == b'\n' || line[n - 1] == b'\r') {
                    n -= 1;
                }
                data.extend_from_slice(&line[..n]);
                stats.bases += n as u64;
                if overlap > 0 && data.len() >= batch_bytes {
                    // cut a long record: carry its last `overlap` bases into the next piece
                    let rec_start = memchr::memrchr(b'\n', &data).map(|i| i + 1).unwrap_or(0);
                    let keep_from = data.len().saturating_sub(overlap).max(rec_start);
                    let tail = data[keep_from..].to_vec();
                    data.push(b'\n');
                    flush!();
                    data.extend_from_slice(&tail);
                }
            }
            if open {
                data.push(b'\n');
            }
        }
        other => {
            return Err(format!("{path}: not FASTQ or FASTA (first byte {:?})", other as char).into());
        }
    }
    flush!();
    Ok(stats)
}

/// Read a whole FASTA of named query sequences (reads or contigs) as ASCII.
pub fn load_named(path: &str) -> Res<Vec<(String, Vec<u8>)>> {
    let mut rdr = open_maybe_gz(path)?;
    let mut out: Vec<(String, Vec<u8>)> = vec![];
    let mut line: Vec<u8> = vec![];
    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            let text = String::from_utf8_lossy(&line[1..]);
            out.push((text.split_whitespace().next().unwrap_or("").to_string(), vec![]));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, content: &[u8]) -> String {
        let p = std::env::temp_dir().join(format!("cpk_test_{}_{}", std::process::id(), name));
        File::create(&p).unwrap().write_all(content).unwrap();
        p.to_string_lossy().to_string()
    }

    #[test]
    fn fasta_packing_roundtrip() {
        let p = tmp("a.fa", b">c1 desc\nACGTNNAC\nGTacgt\n>c2\nNNNN\nA\n");
        let c = load_fasta(&p, |_| Ok(0)).unwrap();
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].name, "c1");
        assert_eq!(c[0].len, 14);
        assert_eq!(c[0].nruns, vec![(4, 6)]);
        let mut codes = vec![];
        c[0].codes(0, 14, &mut codes);
        assert_eq!(codes, vec![0, 1, 2, 3, 4, 4, 0, 1, 2, 3, 0, 1, 2, 3]);
        c[0].codes(5, 9, &mut codes);
        assert_eq!(codes, vec![4, 0, 1, 2]);
        assert_eq!(c[1].len, 5);
        assert_eq!(c[1].nruns, vec![(0, 4)]);
        c[1].codes(0, 5, &mut codes);
        assert_eq!(codes, vec![4, 4, 4, 4, 0]);
        std::fs::remove_file(p).unwrap();
    }

    #[test]
    fn fastq_batches() {
        let p = tmp("r.fq", b"@r1\nACGT\n+\n!!!!\n@r2\nGGNCC\n+\n@@@@@\n");
        let (tx, rx) = crossbeam_channel::unbounded();
        let st = stream_reads(&p, 0, 1 << 20, 0, &tx).unwrap();
        drop(tx);
        let all: Vec<u8> = rx.iter().flat_map(|b| b.data).collect();
        assert_eq!(all, b"ACGT\nGGNCC\n");
        assert_eq!(st.records, 2);
        assert_eq!(st.bases, 9);
        std::fs::remove_file(p).unwrap();
    }

    #[test]
    fn long_fasta_records_are_cut_with_overlap() {
        use crate::kmer::for_each_canon;
        let seq: Vec<u8> = (0..5000u32).map(|i| b"ACGT"[((i * 7 + i / 3 + i * i / 11) % 4) as usize]).collect();
        let mut body = b">chr\n".to_vec();
        for line in seq.chunks(60) {
            body.extend_from_slice(line);
            body.push(b'\n');
        }
        body.extend_from_slice(b">short\nACGTACGTACGTAAA\n");
        let p = tmp("long.fa", &body);
        let k = 11;
        let (tx, rx) = crossbeam_channel::unbounded();
        let st = stream_reads(&p, 0, 700, k - 1, &tx).unwrap();
        drop(tx);
        let batches: Vec<Batch> = rx.iter().collect();
        assert!(batches.len() > 5);
        assert_eq!(st.records, 2);
        assert_eq!(st.bases, 5015);
        let mut got = vec![];
        for b in &batches {
            for_each_canon(&b.data, k, |_, c, _| got.push(c));
        }
        let mut want = vec![];
        for_each_canon(&seq, k, |_, c, _| want.push(c));
        for_each_canon(b"ACGTACGTACGTAAA", k, |_, c, _| want.push(c));
        assert_eq!(got, want); // every k-mer exactly once, in order
        std::fs::remove_file(p).unwrap();
    }

    #[test]
    fn truncated_fastq_is_an_error() {
        let p = tmp("t.fq", b"@r1\nACGT\n+\n!!!!\n@r2\nGG\n");
        let (tx, _rx) = crossbeam_channel::unbounded();
        assert!(stream_reads(&p, 0, 1 << 20, 0, &tx).is_err());
        std::fs::remove_file(p).unwrap();
    }
}
