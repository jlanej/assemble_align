"""Rung 0: a synthetic trio whose assembly carries one planted error of each class."""

import numpy as np
import pandas as pd
import pytest

from copyproof import sim
from copyproof.pipeline import Config, run
from copyproof.reads import ReadSource


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    d = tmp_path_factory.mktemp("case")
    rng = np.random.default_rng(11)
    c = sim.make_case(rng)
    paths = sim.write_case(c, d, rng)
    cfg = Config(
        out=str(d / "out"), asm=[("mat", str(paths["mat"])), ("pat", str(paths["pat"]))],
        reads=ReadSource([str(paths["child"])]),
        parents={"mat": ReadSource([str(paths["mother"])]), "pat": ReadSource([str(paths["father"])])},
        genome_size="4m", threads=4,
    )
    summary = run(cfg)
    issues = pd.read_csv(d / "out" / "issues.tsv", sep="\t")
    return {"dir": d / "out", "summary": summary, "issues": issues, "truth": c["truth"], "case": c, "paths": paths, "cfg": cfg}


def hits(issues, contig, start, end, cls, role="locus"):
    sel = (issues["contig"] == contig) & (issues["class"] == cls) & (issues["role"] == role)
    sel &= (issues["start"] < end) & (issues["end"] > start)
    return issues[sel]


def overlap_fraction(rows, start, end):
    covered = sum(max(0, min(e, end) - max(s, start)) for s, e in zip(rows["start"], rows["end"]))
    return covered / (end - start)


def test_depth_model(case):
    s = case["summary"]
    assert 10 < s["depth_per_copy"] < 13  # 15x per haplotype, 150-bp reads, k = 31
    assert s["assessable_fraction"] > 0.95
    assert s["read_sets"] == ["sample", "mat", "pat"]


def test_collapse_is_found_on_the_right_haplotype(case):
    iss = case["issues"]
    own = hits(iss, "mat", 100_000, 130_000, "COLLAPSE")
    assert len(own) == 1 and overlap_fraction(own, 100_000, 130_000) > 0.95
    assert own.iloc[0]["confidence"] == "high"
    # the homologous stretch is reported as a partner, not as the locus
    assert len(hits(iss, "pat", 100_000, 130_000, "COLLAPSE", role="partner")) == 1
    assert len(hits(iss, "pat", 100_000, 130_000, "COLLAPSE", role="locus")) == 0


def test_missing_junction_of_the_collapsed_duplication(case):
    iss = case["issues"]
    j = iss[iss["class"] == "MISSING_JUNCTION"]
    assert len(j) == 1
    row = j.iloc[0]
    assert row["contig"] == "mat" and abs(row["start"] - 100_000) <= 12
    mate = int(row["evidence"].split("mate=mat:")[1].split(";")[0])
    assert abs(mate - 130_000) <= 12
    # absent from both parents: a de novo junction
    assert "in_parent1=0.0;in_parent2=0.0" in row["evidence"]


def test_expansion_covers_both_copies(case):
    own = hits(case["issues"], "mat", 200_000, 240_000, "EXPANSION")
    assert overlap_fraction(own, 200_000, 240_000) > 0.95
    assert len(hits(case["issues"], "pat", 200_000, 220_000, "EXPANSION", role="partner")) == 1


def test_missing_homologue(case):
    iss = case["issues"]
    own = hits(iss, "mat", 320_000, 335_000, "COLLAPSE")
    assert len(own) == 1 and "holds it once" in own.iloc[0]["summary"]
    hap = hits(iss, "mat", 320_000, 335_000, "MISSING_HAPLOTYPE")
    assert len(hap) == 1 and int(hap.iloc[0]["evidence"].split("alleles=")[1].split(";")[0]) >= 5


def test_haplotype_switch(case):
    iss = case["issues"]
    sw = hits(iss, "mat", 420_000, 450_000, "PHASE_SWITCH")
    assert len(sw) == 1 and overlap_fraction(sw, 420_000, 450_000) > 0.8
    assert len(hits(iss, "mat", 420_000, 450_000, "MISSING_HAPLOTYPE")) == 1
    # a switch moves no copies: the copy-number channel stays silent there
    assert len(iss[(iss["channel"] == "S1") & (iss["contig"] == "mat") & (iss["start"] < 450_000) & (iss["end"] > 420_000)]) == 0


def test_unsupported_sequence(case):
    own = hits(case["issues"], "pat", 485_000, 495_000, "UNSUPPORTED")
    assert len(own) == 1 and overlap_fraction(own, 485_000, 495_000) > 0.95


def test_false_join_where_sequence_was_dropped(case):
    """The assembly's paternal haplotype joins 300 kb to what should follow 15 kb later. No read
    supports that adjacency: the k-mers across it are listed as an unsupported run. (Such runs are
    not issues by themselves; on real data too many correct stretches go unread.)"""
    u = pd.read_csv(case["dir"] / "unsupported_clusters.tsv", sep="\t")
    at = u[(u["contig"] == "pat") & (u["start"] <= 300_000) & (u["end"] >= 300_000)]
    assert len(at) == 1
    row = at.iloc[0]
    assert row["unsupported"] == row["tracked"] >= 4 and row["flank_copies"] > 0.8
    assert row["end"] - row["start"] <= 2 * 31
    assert not (case["issues"]["class"] == "FALSE_JOIN").any()


def test_nothing_else_is_flagged(case):
    iss = case["issues"]
    truth = [(c, s, e) for c, s, e, _ in case["truth"]] + [("pat", 100_000, 130_000), ("pat", 200_000, 220_000), ("pat", 299_000, 301_000)]
    stray = []
    for r in iss.itertuples(index=False):
        if r.role == "twin":
            continue  # the same event, anchored on the other haplotype
        if not any(c == r.contig and r.start < e + 2000 and r.end > s - 2000 for c, s, e in truth):
            stray.append((r.contig, r.start, r.end, r[4]))
    assert stray == []
    # the switched stretch misses the maternal alleles; its flanks also exist on the paternal
    # haplotype, where the same events are twins, not loci. So is the paternal side of the
    # stretch that the paternal haplotype lacks (its reads stop at the false join).
    twins = iss[iss["role"] == "twin"]
    assert len(twins) >= 1 and set(twins["contig"]) == {"pat"}
    assert (twins["start"].between(380_000, 416_000) | twins["start"].between(299_900, 300_100)).all()
    assert (case["dir"] / "issues.bed").read_text().count("\n") == len(iss) + 1


@pytest.fixture(scope="module")
def extra_copy(tmp_path_factory):
    """A correct diploid assembly plus a separate 12-kb sequence that repeats part of the maternal
    haplotype exactly: the commonest kind of false duplication. With no k-mer of its own, the
    extra sequence alone cannot tell "one copy too many" from "two thirds of the depth"."""
    d = tmp_path_factory.mktemp("extra")
    rng = np.random.default_rng(23)
    a, b = sim.diploid(300_000, rng, het=0.001)
    sim.write_fasta(d / "mat.fa", {"mat": a})
    sim.write_fasta(d / "pat.fa", {"pat": b})
    sim.write_fasta(d / "unplaced.fa", {"extra": a[100_000:112_000]})
    sim.write_fastq(d / "child.fq.gz", sim.simulate_reads({"a": a, "b": b}, depth=15, rng=rng, err=0.002))
    cfg = Config(
        out=str(d / "out"), asm=[("mat", str(d / "mat.fa")), ("pat", str(d / "pat.fa")), ("other", str(d / "unplaced.fa"))],
        reads=ReadSource([str(d / "child.fq.gz")]), genome_size="2m", threads=4,
    )
    run(cfg)
    return {"dir": d / "out", "issues": pd.read_csv(d / "out" / "issues.tsv", sep="\t"),
            "segs": pd.read_csv(d / "out" / "cn_segments.tsv", sep="\t")}


def test_identical_extra_copy_is_reported_with_its_copies(extra_copy):
    iss = extra_copy["issues"]
    s1 = iss[iss["channel"] == "S1"]
    assert set(s1["class"]) == {"EXPANSION"}
    ex = s1[s1["contig"] == "extra"]
    assert len(ex) == 1 and ex.iloc[0]["end"] - ex.iloc[0]["start"] >= 10_000
    ev = ex.iloc[0]["evidence"]
    assert "other_copies=" in ev and "mat:100" in ev and "pat:100" in ev
    # the two identical copies cannot be told apart; the homologue is a partner
    mat = s1[(s1["contig"] == "mat") & (s1["start"] < 112_000) & (s1["end"] > 100_000)]
    pat = s1[(s1["contig"] == "pat") & (s1["start"] < 112_000) & (s1["end"] > 100_000)]
    assert len(mat) == 1 and len(pat) == 1 and pat.iloc[0]["role"] == "partner"
    assert ex.iloc[0]["role"] == mat.iloc[0]["role"] == "unresolved"
    assert len(s1) == 3
    # every segment ends where the duplicated stretch ends
    seg = extra_copy["segs"]
    assert len(seg) == 3 and (seg["end"] - seg["start"]).between(11_000, 13_000).all()
    # k-mers unique to the homologue, lying between the shared ones, are read at one copy
    row = seg[seg["contig"] == "extra"].iloc[0]
    n, copies = str(row["copy_unique"]).split(":")
    assert int(n) >= 40 and 0.85 < float(copies) < 1.15 and row["copy_llr"] > 3


def test_copies_settle_an_ambiguous_segment(extra_copy):
    """Had the extra sequence been ambiguous on its own k-mers (as it is on real data, where depth
    is noisier), the unique k-mers of its other copies decide: at full depth it is a surplus
    copy, at the segment's own depth ratio it is a dip in depth."""
    from copyproof.cn import CnParams, assessable
    from copyproof.copies import resolve
    from copyproof.model import fit_model
    from copyproof.windows import Windows

    d = extra_copy["dir"]
    w = Windows.load(str(d / "windows.w1000.npz"))
    model = fit_model(w)
    p = CnParams()
    ok = assessable(w, model, p)
    hits = pd.read_csv(d / "copies.hits.tsv", sep="\t", dtype={"contig": str, "strand": str})
    segs = extra_copy["segs"].copy()
    index = {name: i for i, name in enumerate(w.contigs["name"])}
    segs["w0"] = [int(w.contigs.at[index[c], "win_start"]) + s // w.W for c, s in zip(segs["contig"], segs["start"])]
    segs["w1"] = [int(w.contigs.at[index[c], "win_start"]) + -(-e // w.W) for c, e in zip(segs["contig"], segs["end"])]
    segs["kind"] = "CN_OR_DEPTH"
    segs["llr_depth"] = segs["llr_cn"]  # no evidence of its own either way
    out = resolve(segs, w, model, p, ok, hits, 31)
    row = out[out["contig"] == "extra"].iloc[0]
    assert row["kind"] == "CN" and row["resolved"] == "copies"
    # now deplete the homologue's unique k-mers by the same factor as the segment: a dip in depth
    pat = slice(int(w.contigs.at[index["pat"], "win_start"]) + 100, int(w.contigs.at[index["pat"], "win_start"]) + 112)
    w.S[pat, 1] = np.round(w.S[pat, 1] * 2 / 3).astype(w.S.dtype)
    out = resolve(segs, w, model, p, ok, hits, 31)
    row = out[out["contig"] == "extra"].iloc[0]
    assert row["kind"] == "DEPTH" and row["resolved"] == "copies"


def test_long_reads_turn_unsupported_runs_into_false_joins(case, tmp_path):
    """The same run again with long reads aligned to the assembly. Where 15 kb was dropped, the
    assembly joins sequence that is not adjacent in the genome: no read of either kind holds the
    k-mers across the join and no long read runs through, so the unsupported run becomes a false
    join. The foreign sequence and the false duplication are confirmed: long reads do not carry
    their unread k-mers either."""
    import dataclasses
    import shutil
    import subprocess

    if not (shutil.which("minimap2") and shutil.which("samtools")):
        pytest.skip("needs minimap2 and samtools")
    rng = np.random.default_rng(31)
    truth = case["case"]["child_truth"]
    reads = sim.simulate_reads(truth, depth=20, rng=rng, read_len=10_000, err=0.002)
    (tmp_path / "long.fa").write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
    asm = tmp_path / "asm.fa"
    asm.write_text(open(case["paths"]["mat"]).read() + open(case["paths"]["pat"]).read())
    sam = subprocess.run(["minimap2", "-a", "-x", "map-hifi", "-t", "2", str(asm), str(tmp_path / "long.fa")],
                         capture_output=True, check=True).stdout
    subprocess.run(["samtools", "sort", "-o", str(tmp_path / "long.bam"), "-"], input=sam, check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(tmp_path / "long.bam")], check=True)
    out = tmp_path / "out"
    shutil.copytree(case["dir"], out)
    run(dataclasses.replace(case["cfg"], out=str(out), long_read_bam=str(tmp_path / "long.bam"), cards=0))
    iss = pd.read_csv(out / "issues.tsv", sep="\t")
    fj = iss[iss["class"] == "FALSE_JOIN"]
    assert fj["contig"].tolist() == ["pat"] and abs(int(fj.iloc[0]["start"]) - 300_000) <= 62
    assert fj.iloc[0]["confidence"] == "high" and fj.iloc[0]["channel"] == "S1+S4"
    pts = pd.read_csv(out / "unsupported_points.tsv", sep="\t").set_index("name")
    assert pts.loc["run0", "spanning"] == 0 and pts.loc["run0", "clipped"] >= 5 and pts.loc["run0", "kmers_carried"] == 0
    cn = pts[pts.index.str.startswith("cn")]
    assert len(cn) == 2 and (cn["long_read_kmers"] == "not_carried").all() and (cn["kmers_carried"] == 0).all()
    # the copy-number calls of the first run are all still there, two of them now with long-read support
    before = case["issues"]

    def s1(t):
        t = t[t["channel"].str.startswith("S1") & (t["class"] != "FALSE_JOIN")]
        return sorted(zip(t["contig"], t["class"], t["role"]))

    assert s1(iss) == s1(before)
    confirmed = iss[(iss["channel"] == "S1+S4") & (iss["class"] != "FALSE_JOIN")]
    assert sorted(confirmed["class"]) == ["EXPANSION", "UNSUPPORTED"] and (confirmed["confidence"] == "high").all()
