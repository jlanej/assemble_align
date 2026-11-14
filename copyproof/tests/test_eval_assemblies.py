"""Scoring a run against other assemblies of the same person (scripts/eval_assemblies.py).

The synthetic trio of test_pipeline.py again, scanned with two other "assemblies" as extra read
sets: `good` is the child's true genome; `shared` is the true genome except that it carries the
release's collapse too (the duplicated segment once), as a second assembly can."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from copyproof import sim
from copyproof.pipeline import Config, run
from copyproof.reads import ReadSource

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "eval_assemblies.py"


@pytest.fixture(scope="module")
def scored(tmp_path_factory, cpk):
    d = tmp_path_factory.mktemp("others")
    rng = np.random.default_rng(11)
    c = sim.make_case(rng)
    paths = sim.write_case(c, d, rng)
    child, mother = c["child_truth"], c["mother"]
    sim.write_fasta(d / "good.fa", {"A": child["A"], "B": child["B"]})
    sim.write_fasta(d / "shared.fa", {"A": mother["A"], "B": child["B"]})  # mother's A holds the segment once
    for name in ("good", "shared"):
        subprocess.run([cpk, "count", "--exact", "true", "--out", str(d / name), "--k", "31", "--sample", "4", "--threads", "2",
                        "--genome-size", "4m", str(d / f"{name}.fa")], check=True, capture_output=True)
    cfg = Config(
        out=str(d / "out"), asm=[("mat", str(paths["mat"])), ("pat", str(paths["pat"]))],
        reads=ReadSource([str(paths["child"])]),
        parents={"mat": ReadSource([str(paths["mother"])]), "pat": ReadSource([str(paths["father"])]),
                 "good": str(d / "good.cpkdb"), "shared": str(d / "shared.cpkdb")},
        genome_size="4m", threads=4,
    )
    summary = run(cfg)
    assert summary["read_sets"] == ["sample", "mat", "pat", "good", "shared"]
    # the other assemblies ride along as labels: the calls are those of a run without them
    plain = Config(out=str(d / "plain"), asm=cfg.asm, reads=cfg.reads, genome_size="4m", threads=4,
                   parents={"mat": str(d / "out" / "reads.mat.cpkdb"), "pat": str(d / "out" / "reads.pat.cpkdb")},
                   reads_db=str(d / "out" / "reads.sample.cpkdb"))
    run(plain)
    cols = ["contig", "start", "end", "class", "role", "confidence", "score", "channel", "summary", "evidence"]
    with_labels = pd.read_csv(d / "out" / "issues.tsv", sep="\t")[cols]
    without = pd.read_csv(d / "plain" / "issues.tsv", sep="\t")[cols]
    pd.testing.assert_frame_equal(with_labels, without)
    out = d / "out.eval"
    subprocess.run([sys.executable, str(SCRIPT), "--run", str(d / "out"), "--refs", "good", "shared", "--truth-label", "good",
                    "--out", str(out)], check=True, capture_output=True)
    return {
        "issues": pd.read_csv(f"{out}.cn_issues.tsv", sep="\t"),
        "scores": pd.read_csv(f"{out}.cn_scores.tsv", sep="\t"),
        "segments": pd.read_csv(f"{out}.segments.tsv", sep="\t"),
        "calibration": pd.read_csv(f"{out}.calibration.tsv", sep="\t"),
        "missing": pd.read_csv(f"{out}.missing_clusters.tsv", sep="\t"),
        "unread": pd.read_csv(f"{out}.unread_kmers.tsv", sep="\t").set_index("labels"),
        "json": json.loads(Path(f"{out}.json").read_text()),
    }


def at(issues, contig, start, end):
    return issues[(issues["contig"] == contig) & (issues["start"] < end) & (issues["end"] > start)]


def test_an_error_the_other_assembly_shares_is_confirmed_by_one_only(scored):
    """The collapse: `good` holds one copy more than the release, `shared` holds what the release
    holds. One assembly confirms it, the two together do not."""
    for contig in ("mat", "pat"):  # the locus and its homologous stretch
        c = at(scored["issues"], contig, 100_000, 130_000)
        c = c[c["class"] == "COLLAPSE"]
        assert len(c) == 1
        r = c.iloc[0]
        assert r["confirmed_good"] and not r["confirmed_shared"] and not r["confirmed_both"] and r["confirmed_either"]
        assert r["confirmed_truth"]


def test_errors_of_the_release_alone_are_confirmed_by_both(scored):
    iss = scored["issues"]
    exp = at(iss, "mat", 200_000, 240_000)
    assert len(exp) >= 1 and (exp["class"] == "EXPANSION").all() and exp["confirmed_both"].all()
    junk = at(iss, "pat", 485_000, 495_000)
    assert len(junk) == 1 and junk.iloc[0]["class"] == "UNSUPPORTED" and junk.iloc[0]["confirmed_both"]
    gone = at(iss, "mat", 320_000, 335_000)
    assert len(gone) == 1 and gone.iloc[0]["class"] == "COLLAPSE" and gone.iloc[0]["confirmed_both"]


def test_scores_add_up(scored):
    iss, sc = scored["issues"], scored["scores"]
    row = sc[(sc["direction"] == "both directions") & (sc["size"] == "all")].set_index("labels")
    high = iss[iss["confidence"] == "high"]
    assert row.at["both", "high"] == len(high) == scored["json"]["high_confidence_cn_issues"]
    for name in ("good", "shared", "both", "either", "truth"):
        assert row.at[name, "high_confirmed"] == int(high[f"confirmed_{name}"].sum())
    # every issue is real here, and only the two at the shared collapse escape the stricter rule
    assert row.at["good", "issues_confirmed"] == row.at["good", "issues"] == len(iss)
    assert row.at["both", "issues_confirmed"] == len(iss) - 2
    # the discordances of each label set are found, and "either" has none of its own
    assert row.at["good", "found"] == row.at["good", "assessable"] >= 5
    assert row.at["both", "assessable"] == row.at["good", "assessable"] - 2
    assert pd.isna(row.at["either", "assessable"])
    by_dir = sc[(sc["labels"] == "both") & (sc["size"] == "all")].set_index("direction")
    assert by_dir.at["gain", "issues"] + by_dir.at["loss", "issues"] == by_dir.at["both directions", "issues"]


def test_the_proxy_is_measured_against_the_truth(scored):
    cal = scored["calibration"]
    i = cal[(cal["what"] == "issues") & (cal["group"] == "high")].set_index("labels")
    assert i.at["both", "truth_only"] == 2 and i.at["both", "proxy_only"] == 0
    assert i.at["either", "truth_only"] == 0 and i.at["either", "by_both"] == i.at["either", "n"]
    d = cal[(cal["what"] == "assessable discordances") & (cal["group"] == "all")].set_index("labels")
    assert d.at["both", "proxy_only"] == 0 and d.at["both", "by_truth"] == d.at["both", "n"]


def test_missing_sequence_and_unread_kmers(scored):
    m = scored["missing"]
    assert len(m) and (m["held_by_either"] >= m["held_by_both"]).all() and (m["held_by_good"] <= m["clusters"]).all()
    # the junction of the collapsed duplication is in the true genome and not in the assembly that shares the collapse
    j = m[m["kind"] == "MISSING_JUNCTION"]
    assert int(j["held_by_good"].sum()) == 1 and int(j["held_by_shared"].sum()) == 0
    u = scored["unread"]
    # the 10 kb that is in no genome: unread by the sample and absent from both other assemblies
    assert u.at["both", "unread_and_absent"] >= 2000 and u.at["both", "unread_share_absent"] > 0.9
