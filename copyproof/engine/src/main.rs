//! cpk: the k-mer engine of copyproof.
//!
//!   cpk count    reads -> database of tracked k-mer counts
//!   cpk scan     assembly + read databases -> positional track and absent read k-mers
//!   cpk recruit  reads carrying k-mers of a set
//!   cpk locate   positions of query k-mers in an assembly
//!   cpk extract  sequences of stretches of an assembly
//!   cpk dump     print a database as text

#[macro_use]
mod util;
mod count;
mod fastx;
mod kmer;
mod locate;
mod recruit;
mod scan;
mod table;

const USAGE: &str = "cpk: k-mer engine of copyproof

usage:
  cpk count   --out PREFIX [--k 31] [--sample 4] [--threads N] [--genome-size 3.1g]
              [--min-count 2] [--exact true] [--adapter SEQ|none] READS.fastq[.gz] ... | -
  cpk scan    --out PREFIX --asm LABEL=asm.fa[.gz] [--asm ...] [--group-map names.tsv]
              [--reads LABEL=X.cpkdb ...] [--k 31] [--sample 4] [--novel-min 3] [--threads N]
              [--high-copy 8] [--high-depth 20]
  cpk recruit --kmers KMERS.u64 --out reads.fa [--k 31] [--min-hits 1] [--max-reads N]
              [--repeat-hashes scan.repeat.u64 --sample 4 --max-repeat-frac 0.5] [--set-sample 4]
              READS.fastq[.gz] ... | -
  cpk locate  (--query Q.fa | --regions R.bed) --asm LABEL=asm.fa[.gz] [--asm ...] --out anchors.tsv
  cpk extract --regions R.bed --asm LABEL=asm.fa[.gz] [--asm ...] --out regions.fa
              [--k 31] [--sample 1] [--max-occ 8]
  cpk dump    --db X.cpkdb [--min-count 0]

count cuts reads at the Illumina TruSeq adapter start (AGATCGGAAGAGC) unless --adapter none.
count --exact true counts every sighting (no Bloom filter); use it to count the k-mers of an
assembly, for example a truth genome to compare another assembly with.

A k-mer is tracked when its hash is a multiple of --sample, so every read set and every assembly
tracks the same k-mers. The first --reads set given to scan is the sample itself; any further sets
(parents) are looked up for the same k-mers.";

fn main() {
    util::elapsed();
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 || argv[1] == "-h" || argv[1] == "--help" {
        eprintln!("{USAGE}");
        std::process::exit(if argv.len() < 2 { 2 } else { 0 });
    }
    let a = util::Args::parse(&argv[2..]);
    let r = match argv[1].as_str() {
        "count" => count::run(&a),
        "scan" => scan::run(&a),
        "recruit" => recruit::run(&a),
        "locate" => locate::run(&a),
        "extract" => locate::extract(&a),
        "dump" => count::dump(&a),
        "version" | "--version" => {
            println!("cpk {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        other => Err(format!("unknown command '{other}'\n\n{USAGE}").into()),
    };
    if let Err(e) = r {
        eprintln!("cpk: error: {e}");
        std::process::exit(1);
    }
}
