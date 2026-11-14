"""At short stretches whose haplotype-specific k-mers have no short reads: do HiFi reads carry them?

The first, hand-made version of what copyproof now does with --long-read-bam (NOTEBOOK entry k).
Writes hifi_at_short_losses.tsv."""
import subprocess, sys, os
import pandas as pd
from copyproof.track import Scan
import numpy as np

# hifi_at_short_losses.py OUT_DIR DATA   (DATA: the data root with hg002/ and work/)
D = sys.argv[2]
W = f"{D}/work/hg002_v1.1/run35"
FA = f"{D}/hg002/asm/hg002v1.1.fasta.gz"
BAM = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/hifi_revio_pbmay24/hg002v1.1_hifi_revio_pbmay24.bam"
K = 31
COMP = str.maketrans("ACGT", "TGCA")
def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r

segs = pd.read_csv(f"{W}/cn_segments.tsv", sep="\t")
segs["len"] = segs.end - segs.start
sel = segs[(segs.kind == "CN") & (segs.direction == "loss") & (segs.own == "own") & (segs.len <= 3000) & (segs.contig != "chrM")]
sc = Scan(f"{W}/scan")
names = list(sc.contigs["name"])
os.makedirs(sys.argv[1], exist_ok=True); os.chdir(sys.argv[1])
rows = []
for r in sel.itertuples(index=False):
    ci = names.index(r.contig)
    pos, m, c = sc.region(ci, int(r.start), int(r.end))
    mt = m.sum(axis=1)
    lone = (mt == 1)                       # k-mers found once in the release: this haplotype's own
    unread = lone & (c[:, 0] == 0)
    seq = subprocess.run(["samtools", "faidx", FA, f"{r.contig}:{int(r.start)+1}-{int(r.end)+K-1}"], capture_output=True, text=True, check=True).stdout
    seq = "".join(seq.split("\n")[1:]).upper()
    own = {canon(seq[p - int(r.start): p - int(r.start) + K]): 0 for p in pos[unread]}
    shared = {canon(seq[p - int(r.start): p - int(r.start) + K]): 0 for p in pos[mt == 2]}
    out = subprocess.run(["samtools", "view", "-F", "0x904", BAM, f"{r.contig}:{max(int(r.start)-200,1)}-{int(r.end)+200}"], capture_output=True, text=True)
    reads = [l.split("\t") for l in out.stdout.splitlines()]
    nspan = 0
    for f in reads:
        s = f[9]
        seen = set()
        for i in range(len(s) - K + 1):
            k = canon(s[i:i + K])
            if k in own and k not in seen:
                own[k] += 1; seen.add(k)
            elif k in shared and k not in seen:
                shared[k] += 1; seen.add(k)
    ov = np.array(list(own.values())); sv = np.array(list(shared.values()))
    rows.append({"contig": r.contig, "start": int(r.start), "len": int(r.len), "hifi_reads": len(reads), "own_unread_kmers": len(own),
                 "in_2+_hifi_reads": int((ov >= 2).sum()), "median_hifi_reads": float(np.median(ov)) if len(ov) else float("nan"),
                 "shared_kmers": len(shared), "shared_median_hifi_reads": float(np.median(sv)) if len(sv) else float("nan")})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv("hifi_at_short_losses.tsv", sep="\t", index=False)
