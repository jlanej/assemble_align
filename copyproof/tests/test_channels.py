"""Smaller pieces: intervals, the two-state hapmer decoder, the segment classifier's building
blocks, truth labels, and the long-read join audit on a simulated scaffold."""

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from copyproof import sim
from copyproof.cn import viterbi_star
from copyproof.copies import copy_loci
from copyproof.hapmers import HapParams, audit, pieces_of, readable_windows, viterbi_two_state
from copyproof.hmm import transition_matrix, viterbi
from copyproof.intervals import IntervalSet
from copyproof.missing import low_complexity
from copyproof.kmers import encode
from copyproof.windows import Windows, count_sites, pair_weights, unsupported_clusters


def test_intervals():
    s = IntervalSet.from_rows([("a", 10, 20), ("a", 15, 30), ("a", 50, 60), ("b", 0, 5)])
    assert s.by["a"].tolist() == [[10, 30], [50, 60]] and s.total() == 35 and len(s) == 3
    assert s.overlap("a", 0, 100) == 30 and s.overlap("a", 25, 55) == 10 and s.overlap("c", 0, 9) == 0
    assert np.allclose(s.window_fraction("a", 65, 10), [0, 1, 1, 0, 0, 1, 0])
    u = s.union(IntervalSet.from_rows([("a", 28, 52), ("c", 1, 2)]))
    assert u.by["a"].tolist() == [[10, 60]] and u.total() == 50 + 5 + 1
    assert s.pad(5).by["b"].tolist() == [[0, 10]]


def test_star_viterbi_matches_the_general_one():
    rng = np.random.default_rng(3)
    for _ in range(20):
        n, k = 300, 5
        ll = rng.normal(0, 1, size=(n, k))
        a, b = sorted(rng.integers(0, n, size=2))
        ll[a:b, 2] += 2.5  # a stretch that favours one event state
        p_enter, p_leave = 1e-3, 0.05
        T = transition_matrix(k, p_enter, p_leave)
        start = np.r_[1.0, np.full(k - 1, p_enter / (k - 1))]
        assert viterbi_star(ll, p_enter, p_leave).tolist() == viterbi(ll, T, start).tolist()


def test_two_state_decoder_finds_a_switch():
    a = np.zeros(200, dtype=np.int64)
    b = np.zeros(200, dtype=np.int64)
    a[::5] = 1  # parent 1 sites along the sequence ...
    a[80:120] = 0
    b[80:120:4] = 1  # ... except a stretch that carries parent 2
    b[30] = 1  # one stray site must not make a switch
    path = viterbi_two_state(a, b, err=0.02, p_switch=1e-3)
    assert path[:76].sum() == 0 and path[125:].sum() == 0
    assert path[84:117].min() == 1


def test_pieces_and_sites():
    assert pieces_of(1000, [(100, 200), (500, 505)]) == [(0, 100), (200, 1000)]
    assert pieces_of(1000, [(0, 50), (900, 1000)]) == [(50, 900)]
    pos = np.array([10, 14, 18, 100, 104, 400])
    win = pos // 100
    mask = np.array([True, True, True, True, True, True])
    assert count_sites(pos, win, mask, 31, 5).tolist() == [1, 1, 0, 0, 1]


def test_pair_weights_match_brute_force():
    rng = np.random.default_rng(1)
    key = np.sort(rng.choice(5000, size=400, replace=False)).astype(np.int64)
    brute = np.array([sum(max(0.0, 1 - abs(int(x) - int(y)) / 121) for y in key) for x in key])
    assert np.allclose(pair_weights(key, 121), brute)


def test_unsupported_clusters():
    pos = np.arange(0, 2000, 4)
    zero = np.zeros(len(pos), dtype=bool)
    zero[100:108] = True  # eight tracked k-mers in a row without reads
    zero[300] = True  # a lone one is not a cluster
    copies = np.where(zero, 0.0, 1.0)
    got = unsupported_clusters(pos, zero, copies, 31, 4, "c", mult=np.ones(len(pos), dtype=int))
    assert len(got) == 1
    contig, start, end, n, tracked, flank, nflank, single = got[0]
    assert single == 8
    assert (contig, start, end, n, tracked) == ("c", 400, 428 + 31, 8, 8) and flank == 1.0 and nflank > 50


def test_low_complexity():
    k = 31
    plain = "ACGTTGCAAGGCTTAACGGACGATCGATCGC"
    homo = "ACGTTGCAAGGAAAAAAAACGATCGATCGCA"
    dinuc = "ACGTTGCACACACACACACGATCGATCGCAT"
    assert len(plain) == len(homo) == len(dinuc) == k
    x = np.array([encode(plain), encode(homo), encode(dinuc)], dtype=np.uint64)
    assert low_complexity(x, k).tolist() == [False, True, True]


@pytest.mark.skipif(not (shutil.which("minimap2") and shutil.which("samtools")), reason="needs minimap2 and samtools")
def test_join_audit(tmp_path):
    """Two scaffolds with a 100-N join: one joins truly adjacent sequence, one joins sequence from
    two different molecules. Long reads link the first and not the second."""
    from copyproof.joins import audit_joins, issues_from_joins

    rng = np.random.default_rng(5)
    a, b = sim.random_sequence(60_000, rng), sim.random_sequence(60_000, rng)
    asm = {"good": a[:30_000] + "N" * 100 + a[30_000:], "bad": a[:30_000] + "N" * 100 + b[30_000:]}
    sim.write_fasta(tmp_path / "asm.fa", asm)
    reads = sim.simulate_reads({"a": a, "b": b}, depth=20, rng=rng, read_len=12_000, err=0.002)
    with open(tmp_path / "reads.fa", "w") as f:
        for i, r in enumerate(reads):
            f.write(f">r{i}\n{r}\n")
    sam = subprocess.run(["minimap2", "-a", "-x", "map-hifi", "-t", "2", str(tmp_path / "asm.fa"), str(tmp_path / "reads.fa")],
                         capture_output=True, check=True).stdout
    subprocess.run(["samtools", "sort", "-o", str(tmp_path / "aln.bam"), "-"], input=sam, check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(tmp_path / "aln.bam")], check=True)
    gaps = pd.DataFrame({"contig": ["good", "bad"], "start": [30_000, 30_000], "end": [30_100, 30_100]})
    j = audit_joins(gaps, str(tmp_path / "aln.bam"), tmp_path, lengths={"good": 60_100, "bad": 60_100}).set_index("contig")
    assert j.loc["good", "verdict"] == "linked" and j.loc["good", "linking_reads"] >= 5
    assert j.loc["bad", "verdict"] == "unsupported" and j.loc["bad", "linking_reads"] == 0
    assert j.loc["bad", "clips_left"] + j.loc["bad", "clips_right"] >= 5
    iss = issues_from_joins(j.reset_index())
    assert iss["contig"].tolist() == ["bad"] and iss["class"].tolist() == ["UNSUPPORTED_JOIN"]


def test_copy_loci():
    """Hits of a 10-kb segment at mat:50,000: itself, a copy on the other haplotype with a 40-bp
    deletion, a reversed partial copy, and scattered repeats that are no copy."""
    rng = np.random.default_rng(3)
    q = np.sort(rng.choice(10_000 - 31, 1500, replace=False))
    rows = [("mat", 50_000 + x, "+", x) for x in q]
    rows += [("pat", 70_000 + x - (40 if x > 6_000 else 0), "+", x) for x in q[::2]]
    rows += [("extra", 9_000 - x, "-", x) for x in q[q < 4_000]]
    rows += [("chr9", int(rng.integers(0, 10**6)), "+", x) for x in q[:15]]
    h = pd.DataFrame(rows, columns=["contig", "pos", "strand", "qpos"])
    loci = copy_loci(h, "mat", 50_000, 60_000, 31)
    assert [d["contig"] for d in loci] == ["pat", "extra"]
    pat, extra = loci
    assert pat["strand"] == "+" and abs(pat["start"] - 70_000) < 60 and abs(pat["end"] - 79_960) < 60
    assert pat["kmers"] == len(q[::2])
    assert extra["strand"] == "-" and 5_000 <= extra["start"] < 5_100 and extra["end"] <= 9_031
    assert extra["q_start"] < 100 and 3_900 < extra["q_end"] <= 4_031
    # a tandem copy inside the segment itself is not listed
    assert copy_loci(h[h["contig"] == "mat"], "mat", 50_000, 60_000, 31) == []


@pytest.mark.skipif(not (shutil.which("minimap2") and shutil.which("samtools")), reason="needs minimap2 and samtools")
def test_long_reads_adjudicate_unread_kmers(tmp_path):
    """Assembly k-mers that short reads lack, at three places: correct sequence (long reads carry
    the k-mers), a false join (they do not, and no read runs through), and a wrong base (they do
    not, but reads run through)."""
    from copyproof.joins import adjudicate_cn, audit_points, issues_from_points

    rng = np.random.default_rng(9)
    a, b, c, d = (sim.random_sequence(60_000, rng) for _ in range(4))
    wrong = "A" if d[45_000] != "A" else "C"
    asm = {"right": c, "joined": a[:30_000] + b[30_000:], "typo": d[:45_000] + wrong + d[45_001:]}
    sim.write_fasta(tmp_path / "asm.fa", asm)
    reads = sim.simulate_reads({"a": a, "b": b, "c": c, "d": d}, depth=25, rng=rng, read_len=12_000, err=0.002)
    with open(tmp_path / "reads.fa", "w") as f:
        for i, r in enumerate(reads):
            f.write(f">r{i}\n{r}\n")
    sam = subprocess.run(["minimap2", "-a", "-x", "map-hifi", "-t", "2", str(tmp_path / "asm.fa"), str(tmp_path / "reads.fa")],
                         capture_output=True, check=True).stdout
    subprocess.run(["samtools", "sort", "-o", str(tmp_path / "aln.bam"), "-"], input=sam, check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(tmp_path / "aln.bam")], check=True)
    pts = pd.DataFrame({"name": ["p_right", "p_join", "p_typo"], "contig": ["right", "joined", "typo"],
                        "start": [19_970, 29_970, 44_970], "end": [20_031, 30_031, 45_031]})
    kmers = {n: [asm[c][i : i + 31] for i in range(s + 2, s + 30, 4)] for n, c, s in zip(pts["name"], pts["contig"], pts["start"])}
    out = audit_points(pts, str(tmp_path / "aln.bam"), tmp_path, kmers_of=kmers).set_index("name")
    assert out.loc["p_right", "long_read_kmers"] == "carried" and out.loc["p_right", "kmers_carried"] == 7
    assert out.loc["p_join", "long_read_kmers"] == "not_carried" and out.loc["p_join", "spanning"] == 0
    assert out.loc["p_typo", "long_read_kmers"] == "not_carried" and out.loc["p_typo", "spanning"] >= 3
    iss = issues_from_points(out.reset_index()).set_index("contig")
    assert iss.loc["joined", "class"] == "FALSE_JOIN" and iss.loc["typo", "class"] == "BASE_ERROR" and "right" not in iss.index
    # a short copy-number loss whose k-mers long reads carry is dropped; one they lack becomes firm
    cn = pd.DataFrame({"contig": ["right", "typo"], "start": [19_970, 44_970], "end": [20_031, 45_031], "class": "EXPANSION",
                       "role": "locus", "confidence": "low", "score": 30.0, "channel": "S1", "summary": "x; only 2 assessed windows",
                       "evidence": "delta=-1;low_because=short"})
    kept = adjudicate_cn(cn, out.reset_index())
    assert kept["contig"].tolist() == ["typo"] and kept["confidence"].tolist() == ["high"]
    assert "low_because" not in kept.iloc[0]["evidence"] and "in_long_reads=0" in kept.iloc[0]["evidence"]


def test_hapmer_audit_needs_readable_parents_and_dense_sites():
    """One 1-Mb sequence of the first parent with three stretches that carry the other parent's
    sites: a real switch (dense sites, both parents covered), a stretch the first parent's reads do
    not cover, and a few stray sites spread over 300 kb that has none of its own."""
    nw, cap = 1000, 8
    hap = np.zeros((nw, 3), dtype=np.int32)
    hap[:, 0] = 1  # the sequence's own parent: one site per window
    hap[100:130, 0], hap[100:130, 1] = 0, 2  # a switch
    hap[400:420, 0], hap[400:420, 1] = 0, 2  # the same, where the first parent has no reads
    hap[600:900, 0] = 0
    hap[[620, 700, 780, 860], 1] = 2  # stray sites
    n = np.zeros((nw, cap + 1), dtype=np.int32)
    n[:, 2] = 200
    par = np.full((nw, 2), 200 * 27.0, dtype=np.float32)
    par[400:420, 0] *= 0.1
    contigs = pd.DataFrame({"name": ["seq"], "group": ["hap1"], "length": [nw * 1000], "win_start": [0], "n_win": [nw]})
    zeros = np.zeros((nw, cap + 1), dtype=np.float32)
    w = Windows(W=1000, cap=cap, lam=13.5, contigs=contigs, contig=np.zeros(nw, dtype=np.int32), start=np.arange(nw) * 1000, n=n,
                S=zeros, n_sp=n, S_sp=zeros, zero=n, V=zeros, n_all=n.sum(axis=1), gc=np.full(nw, 0.4, dtype=np.float32),
                acgt=np.full(nw, 1000, dtype=np.int32), hap=hap, par=par)
    gaps = pd.DataFrame(columns=["contig", "start", "end"])
    assert int((~readable_windows(w, HapParams())).sum()) == 20
    _, _, sw = audit(w, gaps, ("mother", "father"))
    assert len(sw) == 1 and (sw.iloc[0]["start"], sw.iloc[0]["end"]) == (100_000, 130_000)
    assert sw.iloc[0]["sites_segment_parent"] == 60 and sw.iloc[0]["segment_parent"] == "father"
    # without the two conditions all three stretches are called
    _, _, raw = audit(w, gaps, ("mother", "father"), HapParams(parent_depth=(0.0, 1e9), min_density=0.0))
    assert sorted(raw["start"]) == [100_000, 400_000, 620_000]
