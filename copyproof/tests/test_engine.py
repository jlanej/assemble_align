"""The engine against a naive pure-Python reference on a small simulated genome."""

import subprocess
from collections import Counter

import numpy as np
import pytest

from copyproof.kmers import canonical_kmers, decode, fmix64, tracked_kmers
from copyproof.track import Scan

K = 31
S = 4


def run(*args, **kw):
    """Run a command, converting every argument to str."""
    return subprocess.run([str(a) for a in args], check=True, capture_output=True, **kw)


def naive_read_counts(reads, k=K, s=S):
    cnt = Counter()
    for r in reads:
        for _, code in tracked_kmers(r, k, s):
            cnt[code] += 1
    return cnt


@pytest.fixture(scope="module")
def counted(cpk, small_diploid):
    d = small_diploid["dir"]
    # A Bloom filter far larger than the data makes false positives negligible, so counts are exact.
    run(cpk, "count", "--out", d / "child", "--k", K, "--sample", S, "--threads", 4, "--adapter", "none",
        "--genome-size", "2m", "--bloom-bits", "400m", d / "r1.fq", d / "r2.fq.gz")
    return d / "child.cpkdb"


def test_count_matches_naive(cpk, small_diploid, counted):
    out = run(cpk, "dump", "--db", counted, text=True).stdout
    got = {}
    for line in out.splitlines():
        kmer, c = line.split("\t")
        got[kmer] = int(c)
    naive = naive_read_counts(small_diploid["reads"])
    want = {decode(code, K): c for code, c in naive.items() if c >= 2}
    assert len(want) > 10_000
    assert got == want


def test_count_reads_stdin_and_agrees(cpk, small_diploid, counted):
    d = small_diploid["dir"]
    with open(d / "r1.fq", "rb") as f1:
        body = f1.read()
    import gzip
    with gzip.open(d / "r2.fq.gz", "rb") as f2:
        body += f2.read()
    run(cpk, "count", "--out", d / "child_stdin", "--k", K, "--sample", S, "--threads", 2, "--adapter", "none",
         "--genome-size", "2m", "--bloom-bits", "400m", "-", input=body)
    a = run(cpk, "dump", "--db", counted, text=True).stdout
    b = run(cpk, "dump", "--db", d / "child_stdin.cpkdb", text=True).stdout
    assert sorted(a.splitlines()) == sorted(b.splitlines())


def test_default_bloom_only_adds_one(cpk, small_diploid):
    """With a small Bloom filter a false positive may add 1 to a count, never more, never less."""
    d = small_diploid["dir"]
    run(cpk, "count", "--out", d / "tight", "--k", K, "--sample", S, "--threads", 4, "--adapter", "none",
         "--genome-size", "2m", "--bloom-bits", "300k", d / "r1.fq", d / "r2.fq.gz")
    out = run(cpk, "dump", "--db", d / "tight.cpkdb", text=True).stdout
    naive = naive_read_counts(small_diploid["reads"])
    truth = {decode(code, K): c for code, c in naive.items()}
    over = 0
    for line in out.splitlines():
        kmer, c = line.split("\t")
        diff = int(c) - truth[kmer]
        assert diff in (0, 1)
        over += diff
    assert over > 0  # the filter is small enough that some false positives occur


def test_scan_track_matches_naive(cpk, small_diploid, counted):
    d = small_diploid["dir"]
    run(cpk, "scan", "--out", d / "scan", "--k", K, "--sample", S, "--threads", 4, "--novel-min", 3,
         "--asm", f"hap1={d / 'hap1.fa'}", "--asm", f"hap2={d / 'hap2.fa.gz'}", "--reads", f"child={counted}")
    sc = Scan(str(d / "scan"))
    assert sc.groups == ["hap1", "hap2"] and sc.reads == ["child"]
    seqs = {**small_diploid["hap1"], **small_diploid["hap2"]}
    group = {"h1_chr": 0, "h2_chr": 1, "h2_extra": 1}
    mult = Counter()
    for name, seq in seqs.items():
        for _, code in tracked_kmers(seq, K, S):
            mult[(code, group[name])] += 1
    naive = naive_read_counts(small_diploid["reads"])
    assert list(sc.contigs["name"]) == ["h1_chr", "h2_chr", "h2_extra"]
    for name, seq in seqs.items():
        pos, m, c = sc.records(name)
        want = list(tracked_kmers(seq, K, S))
        assert pos.tolist() == [p for p, _ in want]
        assert m[:, 0].tolist() == [mult[(code, 0)] for _, code in want]
        assert m[:, 1].tolist() == [mult[(code, 1)] for _, code in want]
        # a read count of 1 is stored as 0: singletons never enter the table
        assert c[:, 0].tolist() == [naive[code] if naive[code] >= 2 else 0 for _, code in want]
        gc, acgt = sc.composition(name)
        assert len(gc) == -(-len(seq) // 100)
        assert int(acgt.sum()) == sum(ch in "ACGT" for ch in seq)
        assert int(gc.sum()) == sum(ch in "GC" for ch in seq)
    # the 100-N gap on haplotype 1 is reported and no k-mer overlaps it
    assert sc.gaps.values.tolist() == [["h1_chr", 15_000, 15_100]]
    pos, _, _ = sc.records("h1_chr")
    assert not np.any((pos > 15_000 - K) & (pos < 15_100))
    # the extra copy of hap2 raises multiplicity there
    pos, m, _ = sc.region("h2_chr", 5_100, 6_900)
    assert (m[:, 1] >= 2).all()


def test_scan_novel_matches_naive(cpk, small_diploid, counted):
    d = small_diploid["dir"]
    sc = Scan(str(d / "scan"))
    kmers, counts = sc.novel()
    seqs = {**small_diploid["hap1"], **small_diploid["hap2"]}
    in_asm = set()
    for seq in seqs.values():
        in_asm.update(code for _, code in tracked_kmers(seq, K, S))
    naive = naive_read_counts(small_diploid["reads"])
    want = {code: c for code, c in naive.items() if c >= 3 and code not in in_asm}
    got = dict(zip(kmers.tolist(), counts[:, 0].tolist()))
    assert got == want
    # the 600 bp missing from both haplotypes are the main source of solid absent k-mers
    a = small_diploid["truth"]["a"]
    hidden = {code for p, code in tracked_kmers(a[30_000 - K + 1 : 30_600 + K - 1], K, S)}
    assert len(hidden) > 100
    assert len(hidden & set(got)) > 0.9 * len(hidden)
    solid = {code for code, c in got.items() if c >= 10}
    assert len(solid - hidden) < 0.1 * len(solid)


def test_recruit_and_locate(cpk, small_diploid, counted, tmp_path):
    d = small_diploid["dir"]
    sc = Scan(str(d / "scan"))
    kmers, counts = sc.novel()
    solid = kmers[counts[:, 0] >= 8]
    assert len(solid) > 5
    solid.astype("<u8").tofile(tmp_path / "set.u64")
    run(cpk, "recruit", "--kmers", tmp_path / "set.u64", "--out", tmp_path / "rec.fa", "--k", K, "--adapter", "none",
         "--threads", 2, d / "r1.fq", d / "r2.fq.gz")
    solid_set = set(solid.tolist())
    want = 0
    for r in small_diploid["reads"]:
        if any(code in solid_set for _, code, _ in canonical_kmers(r, K)):
            want += 1
    lines = (tmp_path / "rec.fa").read_text().splitlines()
    names, seqs = lines[0::2], lines[1::2]
    assert len(names) == want > 0
    # every listed hit is a k-mer of the set at the stated read position
    for name, seq in zip(names[:50], seqs[:50]):
        hits = name.split("k=")[1].split(",")
        for h in hits:
            ki, p = map(int, h.split(":"))
            sub = seq[p : p + K]
            from copyproof.kmers import canonical
            assert canonical(sub) == int(solid[ki])
    # locate: recruited reads reach into the flanks of the missing 600 bp, which the assembly holds
    run(cpk, "locate", "--query", tmp_path / "rec.fa", "--out", tmp_path / "anchors.tsv", "--k", K,
         "--asm", f"hap1={d / 'hap1.fa'}", "--asm", f"hap2={d / 'hap2.fa.gz'}")
    import pandas as pd
    an = pd.read_csv(tmp_path / "anchors.tsv", sep="\t")
    assert set(an["contig"]) <= {"h1_chr", "h2_chr", "h2_extra"}
    b = small_diploid["hap2"]["h2_chr"]
    assert ((an["contig"] == "h2_chr") & (an["pos"].between(29_800, 30_200))).sum() > 50
    row = an[an["contig"] == "h2_chr"].iloc[0]
    q = dict(zip([n[1:].split()[0] for n in names], seqs))[row["query"]]
    kq = q[row["qpos"] : row["qpos"] + K]
    ka = b[row["pos"] : row["pos"] + K]
    from copyproof.sim import revcomp
    assert kq == (ka if row["strand"] == "+" else revcomp(ka))


def test_hash_matches_engine():
    # the same constants as engine/src/kmer.rs; sampling must agree bit for bit
    assert fmix64(0) == 0
    assert fmix64(1) == 0xB456BCFC34C2CB2C


def test_exact_count_of_an_assembly(cpk, small_diploid, tmp_path):
    """--exact counts every sighting, so an assembly's k-mer multiplicities come out exactly,
    including for sequences far longer than one batch."""
    d = small_diploid["dir"]
    run(cpk, "count", "--out", tmp_path / "asm", "--k", K, "--sample", S, "--threads", 3, "--exact", "true",
        "--genome-size", "2m", d / "hap1.fa", d / "hap2.fa.gz")
    out = run(cpk, "dump", "--db", tmp_path / "asm.cpkdb", text=True).stdout
    got = {}
    for line in out.splitlines():
        kmer, c = line.split("\t")
        got[kmer] = int(c)
    want = Counter()
    for seq in list(small_diploid["hap1"].values()) + list(small_diploid["hap2"].values()):
        for _, code in tracked_kmers(seq, K, S):
            want[decode(code, K)] += 1
    assert got == dict(want)
    assert min(got.values()) == 1 and max(got.values()) >= 3


def test_adapters_are_trimmed(cpk, small_diploid, tmp_path):
    """Reads from fragments shorter than the read run into the sequencing adapter. With trimming
    (the default) the counts equal those of the inserts alone."""
    import random

    rnd = random.Random(5)
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    inserts, reads = [], []
    for r in small_diploid["reads"][:4000]:
        n = rnd.choice([60, 90, 120, 140, 150])
        # skip the rare insert that itself ends like the adapter's start: trimming would shorten it
        ins = r[:n]
        if any(ins.endswith(adapter[:j]) for j in range(5, 14)) or adapter[:13] in ins:
            continue
        inserts.append(ins)
        reads.append((ins + adapter)[:150])
    with open(tmp_path / "ad.fq", "w") as f:
        for i, r in enumerate(reads):
            f.write(f"@a{i}\n{r}\n+\n{'I' * len(r)}\n")
    naive = naive_read_counts(inserts)
    want = {decode(code, K): c for code, c in naive.items() if c >= 2}
    run(cpk, "count", "--out", tmp_path / "trim", "--k", K, "--sample", S, "--threads", 2,
        "--genome-size", "2m", "--bloom-bits", "400m", tmp_path / "ad.fq")
    out = run(cpk, "dump", "--db", tmp_path / "trim.cpkdb", text=True).stdout
    got = {line.split("\t")[0]: int(line.split("\t")[1]) for line in out.splitlines()}
    assert got == want and len(want) > 300
    # without trimming, adapter k-mers appear with high counts
    run(cpk, "count", "--out", tmp_path / "raw", "--k", K, "--sample", S, "--threads", 2, "--adapter", "none",
        "--genome-size", "2m", "--bloom-bits", "400m", tmp_path / "ad.fq")
    out = run(cpk, "dump", "--db", tmp_path / "raw.cpkdb", text=True).stdout
    raw = {line.split("\t")[0]: int(line.split("\t")[1]) for line in out.splitlines()}
    assert max(raw.values()) > 500 and max(got.values()) < 100


def test_locate_regions_lists_the_other_copies(cpk, small_diploid, tmp_path):
    """A stretch of the release as the query: the 2-kb extra sequence is a copy of h2_chr
    5,000-7,000, and shares most of its k-mers with the other haplotype there."""
    import pandas as pd

    d = small_diploid["dir"]
    (tmp_path / "r.bed").write_text("h2_extra\t0\t2000\tseg7\n")
    run(cpk, "locate", "--regions", tmp_path / "r.bed", "--out", tmp_path / "hits.tsv", "--k", K, "--sample", S,
        "--asm", f"hap1={d / 'hap1.fa'}", "--asm", f"hap2={d / 'hap2.fa.gz'}")
    h = pd.read_csv(tmp_path / "hits.tsv", sep="\t")
    assert set(h["query"]) == {"seg7"}
    extra = small_diploid["hap2"]["h2_extra"]
    want = {p for p, _ in tracked_kmers(extra, K, S)}
    own = h[h["contig"] == "h2_extra"]
    assert set(own["qpos"]) == want and (own["qpos"] == own["pos"]).all() and (own["strand"] == "+").all()
    twin = h[h["contig"] == "h2_chr"]
    assert set(twin["qpos"]) == want and (twin["pos"] - twin["qpos"] == 5000).all()
    hom = h[h["contig"] == "h1_chr"]
    assert 0.5 * len(want) < len(hom) < len(want)  # the other haplotype differs at heterozygous sites
    assert (hom["pos"] - hom["qpos"] == 5000).all()
    assert (h.groupby("qpos").size() == h.groupby("qpos")["occ"].first()).all()
