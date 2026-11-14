"""Run every stage for one sample and write its outputs."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from . import calls as C
from .cn import CnParams, call_segments
from .copies import locate_segments
from .copies import resolve as resolve_with_copies
from .engine import engine_path
from .hapmers import HapParams, audit
from .missing import MissingParams, find_missing
from .model import Model, fit_model
from .reads import ReadSource
from .track import Scan
from .windows import UNSUPPORTED_COLUMNS, Windows, compute_windows


@dataclass
class Config:
    out: str
    asm: list[tuple[str, str]]  # (group label, FASTA)
    reads: ReadSource | None = None  # the sample's reads ...
    reads_db: str | None = None  # ... or an existing database of them
    parents: dict = field(default_factory=dict)  # name -> ReadSource or path of a database
    group_regex: list[tuple[str, str]] = field(default_factory=list)  # (label, regex on sequence names)
    group_map: str | None = None
    k: int = 31
    sample: int = 4
    threads: int = 8
    window: int = 1000
    cap: int = 8
    genome_size: str = "3.1g"
    missing: bool = True
    cards: int = 200  # evidence cards to draw, best-scoring first (0: none)
    long_read_bam: str | None = None  # long reads aligned to the release (path or URL, indexed): join audit
    max_points: int = 3000  # unsupported runs checked against the long reads, strongest first
    cn: CnParams = field(default_factory=CnParams)
    hap: HapParams = field(default_factory=HapParams)
    miss: MissingParams = field(default_factory=MissingParams)


def _log(msg: str) -> None:
    print(f"[copyproof {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def fasta_names(path: str) -> list[str]:
    """Sequence names of a FASTA, from its .fai when present."""
    fai = Path(str(path) + ".fai")
    if fai.is_file():
        return [line.split("\t")[0] for line in fai.read_text().splitlines() if line]
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def build_group_map(cfg: Config, path: Path) -> str | None:
    """Assign sequences to groups by regular expression (first match wins)."""
    if cfg.group_map:
        return cfg.group_map
    if not cfg.group_regex:
        return None
    rules = [(label, re.compile(rx)) for label, rx in cfg.group_regex]
    with open(path, "w") as out:
        for _, fa in cfg.asm:
            for name in fasta_names(fa):
                for label, rx in rules:
                    if rx.search(name):
                        out.write(f"{name}\t{label}\n")
                        break
    return str(path)


def _count(source: ReadSource, prefix: Path, cfg: Config, log) -> str:
    db = Path(f"{prefix}.cpkdb")
    if not db.is_file():
        _log(f"counting k-mers -> {db.name}")
        tmp = Path(f"{prefix}.tmp")
        source.run(["count", "--out", tmp, "--k", cfg.k, "--sample", cfg.sample, "--threads", cfg.threads,
                    "--genome-size", cfg.genome_size], log=log)
        for ext in (".cpkdb", ".hist.tsv", ".json"):
            Path(f"{tmp}{ext}").rename(f"{prefix}{ext}")
    return str(db)


def _md5_head(path: str, nbytes: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read(nbytes))
    return h.hexdigest()


def run(cfg: Config) -> dict:
    """Run copyproof on one sample. Stages whose outputs exist are reused."""
    t0 = time.time()
    out = Path(cfg.out)
    out.mkdir(parents=True, exist_ok=True)
    log = open(out / "engine.log", "a")
    timings = {}

    # 1. k-mer databases
    t = time.time()
    labels, dbs = [], []
    if cfg.reads_db:
        child_db = cfg.reads_db
    elif cfg.reads is not None:
        child_db = _count(cfg.reads, out / "reads.sample", cfg, log)
    else:
        raise ValueError("give the sample's reads or a database of them")
    labels.append("sample")
    dbs.append(child_db)
    for name, src in cfg.parents.items():
        db = src if isinstance(src, str) else _count(src, out / f"reads.{name}", cfg, log)
        labels.append(name)
        dbs.append(db)
    timings["count_s"] = round(time.time() - t, 1)

    # 2. scan
    t = time.time()
    prefix = out / "scan"
    asm_args = []
    for label, fa in cfg.asm:
        asm_args += ["--asm", f"{label}={fa}"]
    gm = build_group_map(cfg, out / "groups.tsv")
    if gm:
        asm_args += ["--group-map", gm]
    if not Path(f"{prefix}.scan.json").is_file():
        _log("scanning the release")
        args = [engine_path(), "scan", "--out", str(prefix), "--k", str(cfg.k), "--sample", str(cfg.sample),
                "--threads", str(cfg.threads), "--high-copy", str(cfg.cap)] + asm_args
        for label, db in zip(labels, dbs):
            args += ["--reads", f"{label}={db}"]
        subprocess.run(args, check=True, stderr=log)
    sc = Scan(str(prefix))
    timings["scan_s"] = round(time.time() - t, 1)

    # 3. windows and the depth model
    t = time.time()
    parents = (1, 2) if len(sc.reads) >= 3 else None
    wpath = out / f"windows.w{cfg.window}.npz"
    upath = out / "unsupported_clusters.tsv"
    if wpath.is_file() and upath.is_file():
        w = Windows.load(str(wpath))
        unsupported = pd.read_csv(upath, sep="\t")
    else:
        _log("aggregating windows")
        clusters = []
        w = compute_windows(sc, W=cfg.window, cap=cfg.cap, parents=parents, unsupported=clusters,
                            processes=max(1, min(cfg.threads, 6)))
        w.save(str(wpath))
        unsupported = pd.DataFrame(clusters, columns=UNSUPPORTED_COLUMNS)
        unsupported.to_csv(upath, sep="\t", index=False)
    model = fit_model(w)
    (out / "model.json").write_text(model.to_json())
    timings["windows_s"] = round(time.time() - t, 1)

    # 4. S1 copy number
    t = time.time()
    _log("S1: copy number")
    segs, ok = call_segments(w, model, cfg.cn)
    # The other copies of each candidate segment name the sequences involved, and the k-mers
    # unique to them settle segments that a depth factor explains as well as a whole copy.
    cand = segs[segs["kind"].isin(["CN", "CN_OR_DEPTH"])]
    if len(cand):
        hits = locate_segments(cand, asm_args, sc.k, sc.sample, cfg.cap, out / "copies", cfg.threads, log)
        segs = resolve_with_copies(segs, w, model, cfg.cn, ok, hits, sc.k)
    segs.drop(columns=["ci", "w0", "w1"]).to_csv(out / "cn_segments.tsv", sep="\t", index=False)
    C.depth_anomalies(segs).to_csv(out / "depth_anomalies.tsv", sep="\t", index=False)
    mask = C.mask_runs(w, ok)
    mask.to_csv(out / "cannot_assess.bed", sep="\t", index=False, header=False)
    single_group = sum(1 for g in sc.meta["groups"] if g["sequences"] > 0) < 2
    tables = []
    timings["cn_s"] = round(time.time() - t, 1)

    # 5. S3 hapmers
    hap_summary = {}
    if parents is not None:
        t = time.time()
        _log("S3: parental hapmers")
        names = (sc.reads[1], sc.reads[2])
        contigs, pieces, switches = audit(w, sc.gaps, names, cfg.hap)
        contigs.to_csv(out / "hapmer_sequences.tsv", sep="\t", index=False)
        pieces.to_csv(out / "hapmer_pieces.tsv", sep="\t", index=False)
        switches.to_csv(out / "hapmer_switches.tsv", sep="\t", index=False)
        tables.append(C.issues_from_hapmers(switches, pieces))
        a, b = contigs[f"sites_{names[0]}"].sum(), contigs[f"sites_{names[1]}"].sum()
        assessed = contigs[contigs["parent"] != "not_assessable"]
        hap_summary = {
            "hapmer_sites": int(a + b),
            "sequences_assessed": int(len(assessed)),
            "sequences_not_assessable": int(len(contigs) - len(assessed)),
            "switch_segments": int(len(switches)),
            "switched_bp": int(assessed["switched_bp"].sum()),
            "pieces_other_parent": int((pieces["call"] == "other_parent").sum()) if len(pieces) else 0,
        }
        timings["hapmers_s"] = round(time.time() - t, 1)

    # 6. S2 missing k-mers
    miss_summary = {}
    clusters = pd.DataFrame()
    if cfg.missing and cfg.reads is not None:
        t = time.time()
        _log("S2: read k-mers missing from the release")
        miss_summary, clusters, reads = find_missing(sc, model.lam, cfg.reads, asm_args, out / "missing", cfg.miss, parents,
                                                     cfg.threads, log)
        if len(clusters):
            clusters.to_csv(out / "missing_clusters.tsv", sep="\t", index=False)
            # where both haplotypes could host a missing allele, prefer the one another channel points to
            firm = C.issues_from_cn(segs, single_group)
            hints = [t for t in tables + [firm[firm["role"] == "locus"]] if len(t)]
            tables.append(C.issues_from_missing(clusters, sc.sample, prefer=pd.concat(hints, ignore_index=True) if hints else None))
        timings["missing_s"] = round(time.time() - t, 1)
    else:
        from .missing import solid_missing, summarise
        kmers, counts, thr, cand = solid_missing(sc, model.lam, cfg.miss)
        miss_summary = summarise(sc, model.lam, kmers, counts, thr, cand, parents, cfg.miss)
        miss_summary["note"] = "missing k-mers were counted but not localised (no reads given, or --no-missing)"

    # 7. copy-number issues, now that the other channels can vouch for the ambiguous segments
    corroborated = C.corroborate(segs, tables)
    cn_issues = C.issues_from_cn(segs, single_group, corroborated, cfg.cn.deep_ratio, cfg.cn.deep_windows,
                                 cfg.cn.high_windows)
    tables.insert(0, cn_issues)

    # 8. assembly k-mers without read support: wrong bases and false joins
    alleles = clusters if (cfg.missing and cfg.reads is not None and len(clusters)) else None
    tables.append(C.issues_from_unsupported(unsupported, alleles, sc.k, w, ok, covered=cn_issues))

    # 9. S4: does any long read link the two sides of each run of N?
    join_summary = {}
    if cfg.long_read_bam:
        from .joins import audit_joins, issues_from_joins

        t = time.time()
        _log("S4: join audit")
        lengths = dict(zip(sc.contigs["name"], sc.contigs["length"]))
        joins = audit_joins(sc.gaps, cfg.long_read_bam, out / "joins", lengths=lengths)
        joins.to_csv(out / "joins.tsv", sep="\t", index=False)
        join_at = len(tables)
        tables.append(issues_from_joins(joins))
        # an ambiguous copy-number segment right at an unlinked join is worth reporting
        at_join = C.segments_at_joins(segs, joins)
        if at_join:
            cn_issues = C.issues_from_cn(segs, single_group, corroborated, cfg.cn.deep_ratio, cfg.cn.deep_windows,
                                         cfg.cn.high_windows, at_join=at_join)
            tables[0] = cn_issues
        join_summary = {k: int(v) for k, v in joins["verdict"].value_counts().items()} if len(joins) else {}
        # Assembly k-mers that no short read holds: runs with normally covered flanks, and
        # copy-number losses that sit on this haplotype's own k-mers. Do long reads hold them?
        from .joins import adjudicate_cn, audit_points, extract_regions, issues_from_points, unread_kmers

        u = unsupported[(unsupported["unsupported"] >= unsupported["tracked"]) & (unsupported["flank_kmers"] >= 10)
                        & (unsupported["flank_copies"] >= 0.85)]
        u = u.sort_values("unsupported", ascending=False).head(cfg.max_points)
        runs = pd.DataFrame({"name": [f"run{i}" for i in range(len(u))], "contig": u["contig"].values,
                             "start": u["start"].values, "end": u["end"].values + sc.k - 1})
        lost = cn_issues[(cn_issues["class"].isin(["EXPANSION", "UNSUPPORTED"])) & (cn_issues["role"] == "locus")]
        own = pd.DataFrame({"name": [f"cn{i}" for i in range(len(lost))], "contig": lost["contig"].values,
                            "start": lost["start"].values, "end": lost["end"].values})
        targets = pd.concat([runs, own], ignore_index=True)
        if len(targets):
            _log(f"S4: long reads at {len(runs)} runs of unsupported assembly k-mers and {len(own)} copy-number losses")
            # thousands of remote fetches: reuse them when the same places are asked of the same file
            stamp = hashlib.md5((cfg.long_read_bam + targets.to_csv(index=False)).encode()).hexdigest()
            stamp_file, ppath = out / "joins" / "points.stamp", out / "unsupported_points.tsv"
            points = None
            if ppath.is_file():
                old = pd.read_csv(ppath, sep="\t", dtype={"name": str, "contig": str})
                cols = ["name", "contig", "start", "end"]
                same = (len(old) == len(targets) and "long_read_kmers" in old.columns
                        and (old[cols].astype(str).values == targets[cols].astype(str).values).all())
                if same and (not stamp_file.is_file() or stamp_file.read_text() == stamp):
                    points = old
            if points is None:
                seqs = extract_regions(targets.assign(end=targets["end"] + sc.k), asm_args, out / "joins" / "points", log)
                kmers = unread_kmers(sc, targets, seqs, cfg.cap)
                points = audit_points(targets, cfg.long_read_bam, out / "joins", threads=cfg.threads, kmers_of=kmers)
                points.to_csv(ppath, sep="\t", index=False)
                stamp_file.write_text(stamp)
            is_run = points["name"].str.startswith("run")
            tables.append(issues_from_points(points[is_run]))
            tables[0] = adjudicate_cn(tables[0], points[~is_run])
            join_summary["unsupported_runs"] = {k: int(v) for k, v in points.loc[is_run, "long_read_kmers"].value_counts().items()}
            join_summary["copy_number_losses"] = {k: int(v) for k, v in points.loc[~is_run, "long_read_kmers"].value_counts().items()}
        # an unlinked join is firm only next to another issue
        from .joins import corroborate_joins

        tables[join_at] = corroborate_joins(tables[join_at], [t for i, t in enumerate(tables) if i != join_at])
        timings["joins_s"] = round(time.time() - t, 1)

    # 10. issues and summary
    issues = C.combine(*tables)
    issues.to_csv(out / "issues.tsv", sep="\t", index=False)
    C.write_bed(issues, out / "issues.bed")
    if cfg.cards > 0 and len(issues):
        from .cards import write_cards

        # single alleles and wrong bases are too many and too small for a card each
        pick = issues[~issues["class"].isin(["MISSING_ALLELE", "BASE_ERROR"])]
        pick = pick.sort_values(["confidence", "score"], ascending=[True, False]).head(cfg.cards)
        _log(f"drawing {len(pick)} evidence cards")
        write_cards(issues, w, model, ok, sc.gaps, out / "cards", (sc.reads[1], sc.reads[2]) if parents else None, only=pick)
    total = int(w.contigs["length"].sum())
    assessed_bp = int(ok.sum()) * w.W
    sp = sc.spectrum()
    solid_in = float(sp.loc[sp["count"] >= miss_summary.get("threshold", cfg.miss.min_count), "n_kmers"].sum()) if miss_summary else 0.0
    solid_out = float(miss_summary.get("tracked_kmers", 0))
    by_class = issues.groupby(["class", "role", "confidence"]).size().reset_index(name="n")
    summary = {
        "copyproof": __version__,
        "k": sc.k, "sample": sc.sample, "window": w.W, "multiplicity_cap": w.cap,
        "groups": sc.meta["groups"],
        "read_sets": sc.reads,
        "depth_per_copy": round(model.lam, 3),
        "depth_model": {"omega": model.omega, "psi": round(model.psi, 4), **model.notes},
        "bases": total,
        "assessable_fraction": round(assessed_bp / max(total, 1), 4),
        "cannot_assess_bp": int(mask["end"].sub(mask["start"]).sum()) if len(mask) else 0,
        "kmer_support": C.kmer_support(sp, sc.k, sc.sample),
        "read_kmer_completeness": round(solid_in / (solid_in + solid_out), 6) if solid_in + solid_out else None,
        "cn_segments": {k: int(v) for k, v in segs["kind"].value_counts().items()},
        "cn_or_depth_not_reported": int((segs["kind"] == "CN_OR_DEPTH").sum()
                                        - cn_issues["evidence"].str.contains("low_because=depth_factor_fits").sum()),
        "cn_settled_by_other_copies": int((segs["resolved"] == "copies").sum()) if "resolved" in segs.columns else 0,
        "unsupported_clusters": int(len(unsupported)),
        "hapmers": hap_summary,
        "missing": miss_summary,
        "joins": join_summary,
        "issues": int(len(issues)),
        "issues_by_class": by_class.to_dict("records"),
        "timings": {**timings, "total_s": round(time.time() - t0, 1)},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    run_info = {
        "copyproof": __version__,
        "engine": subprocess.run([engine_path(), "version"], capture_output=True, text=True).stdout.strip(),
        "asm": [{"label": l, "path": str(p), "md5_first_mb": _md5_head(str(p))} for l, p in cfg.asm],
        "databases": dict(zip(labels, dbs)),
        "parameters": {"k": cfg.k, "sample": cfg.sample, "window": cfg.window, "cap": cfg.cap, "cn": asdict(cfg.cn),
                       "hapmers": asdict(cfg.hap), "missing": asdict(cfg.miss)},
        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (out / "run.json").write_text(json.dumps(run_info, indent=1, default=str))
    _log(f"done: {len(issues)} issues; {summary['assessable_fraction'] * 100:.1f}% of bases assessable")
    return summary
