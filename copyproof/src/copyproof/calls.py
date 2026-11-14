"""Turn the evidence channels into one list of issues, a mask and summary numbers.

An issue is a flagged locus with its evidence. It becomes an error only after adjudication.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .windows import Windows

# BED itemRgb per class
COLOURS = {
    "COLLAPSE": "214,96,0",
    "EXPANSION": "0,114,178",
    "UNSUPPORTED": "86,180,233",
    "PHASE_SWITCH": "204,121,167",
    "MISPLACED_PIECE": "170,51,119",
    "MISSING_JUNCTION": "213,94,0",
    "MISSING_ALLELE": "230,159,0",
    "MISSING_HAPLOTYPE": "0,158,115",
    "MISSING_SEQUENCE": "240,228,66",
    "UNPLACED": "128,128,128",
    "FALSE_JOIN": "220,5,12",
    "UNSUPPORTED_JOIN": "136,34,85",
    "BASE_ERROR": "153,153,153",
}

ISSUE_COLUMNS = [
    "id", "contig", "start", "end", "class", "role", "confidence", "score", "channel", "summary", "evidence",
]


# Sequences whose copy number per cell is not one per haplotype: mitochondrion, EBV episome.
EPISOME = r"(?i)(^|[_#|:.-])(chr)?(m|mt|ebv)$|mito|ebv"


def long_and_deep(r, deep_ratio: float = 1.5, deep_windows: int = 8) -> bool:
    """Does the whole-copy reading of a segment change the depth of its commonest class by at
    least `deep_ratio`, over at least `deep_windows` assessed windows? On a truth genome almost
    no ambiguous segment does (short dips to half the depth are common; long ones are not)."""
    m, d = float(r.dominant_m), float(r.delta)
    if m < 1 or d == 0 or int(r.n_assessed) < deep_windows:
        return False
    after = max(m + d, 0.0)
    return after * deep_ratio <= m + 1e-9 or after >= m * deep_ratio - 1e-9


def issues_from_cn(segs: pd.DataFrame, single_group: bool = False, corroborated: set | None = None,
                   deep_ratio: float = 1.5, deep_windows: int = 8, high_windows: int = 5,
                   at_join: set | None = None) -> pd.DataFrame:
    """Candidate issues from the copy-number channel. Segments explained by a depth factor alone,
    or by nothing regular, are not issues (see depth_anomalies). With a single assembly group
    (haplotypes not labelled) there are no haplotype-specific k-mers, so every issue's role is
    `unresolved`: the change is here or on another copy of the same sequence.

    A CN_OR_DEPTH segment (a depth factor fits it as well as a whole copy) is an issue, of low
    confidence, only with support: k-mers of its own haplotype that follow the change, another
    channel (`corroborated` holds the index of those, see `corroborate`), or length and depth
    (`long_and_deep`). The rest are left out, because on a truth genome nearly all of them are
    short dips of sequencing depth. A whole-copy call on fewer than `high_windows` assessed
    windows is of low confidence too (`low_because=short`)."""
    import re

    rows = []
    epi = re.compile(EPISOME)
    has_copies = "copies" in segs.columns
    for i, r in zip(segs.index, segs.itertuples(index=False)):
        if r.kind not in ("CN", "CN_OR_DEPTH", "ABSENT"):
            continue
        if epi.search(str(r.contig)):
            continue  # many copies per cell by nature
        if r.kind == "ABSENT":
            rows.append(
                {
                    "contig": r.contig, "start": int(r.start), "end": int(r.end), "class": "UNSUPPORTED",
                    "role": "locus", "confidence": "high", "score": float(r.llr), "channel": "S1",
                    "summary": f"the reads hold at most {r.factor:g} of the expected depth for every k-mer class here: "
                               "sequence that most or all of the reads' DNA lacks",
                    "evidence": f"factor={r.factor:g};copies_by_multiplicity={r.classes};fit_depth={r.fit_depth};kmers={r.n_kmers};"
                                f"windows={r.n_assessed}",
                }
            )
            continue
        # A segment that a depth factor explains as well as a copy change is an issue only with
        # support of its own: haplotype-specific k-mers that follow the change, or another channel.
        why = ""
        if r.kind == "CN_OR_DEPTH":
            if r.own == "own":
                why = "own_kmers"
            elif corroborated is not None and i in corroborated:
                why = "other_channel"
            elif long_and_deep(r, deep_ratio, deep_windows):
                why = "long_and_deep"
            elif at_join is not None and i in at_join:
                why = "unlinked_join"  # a copy too many or too few right where no long read links the scaffold
            else:
                continue
        gain = r.delta > 0
        only_copy = r.dominant_m == 1
        if gain:
            cls = "COLLAPSE"
            what = (
                f"reads hold {abs(r.delta)} more cop{'y' if abs(r.delta) == 1 else 'ies'} of this sequence than the release"
            )
            if only_copy and r.own == "own":
                what += "; the release holds it once (a collapsed duplicate, or the other haplotype's copy is missing)"
        elif only_copy and r.own == "own":
            cls = "UNSUPPORTED"
            what = "sequence found once in the release has no read support"
        else:
            cls = "EXPANSION"
            what = "the release holds one more copy of this sequence than the reads support"
        role = {"own": "locus", "other": "partner", "unresolved": "unresolved", "na": "unresolved"}.get(r.own, "unresolved")
        if single_group:
            role = "unresolved"
        if role == "partner":
            what += "; k-mers unique to this haplotype are concordant, so the change is on the homologous or a paralogous copy"
        short = r.kind == "CN" and int(r.n_assessed) < high_windows
        conf = "high" if r.kind == "CN" and not short else "low"
        copies = str(r.copies) if has_copies and isinstance(r.copies, str) else ""
        resolved = has_copies and r.resolved == "copies"
        if r.kind == "CN_OR_DEPTH":
            what += "; a depth factor of " + f"{r.factor:g} explains it equally well (a deletion or gain in the reads' DNA, or uneven depth)"
        elif resolved:
            what += "; k-mers unique to the other copies of this sequence are read at full depth, which a dip in sequencing depth would not leave"
        if short:
            what += f"; only {int(r.n_assessed)} assessed windows, and short stretches that short reads fail on look the same"
        if copies and role == "unresolved" and not gain:
            first = copies.split(",")[0].rsplit(":", 1)[0]
            what += f"; the surplus copy is this one or one it shares its k-mers with ({first})"
        ev = (f"delta={r.delta};copies_by_multiplicity={r.classes};own_specific={getattr(r, 'own_classes', '')};"
              f"own_llr={r.own_llr};fit_cn={r.fit_cn};fit_depth={r.fit_depth};kmers={r.n_kmers};windows={r.n_assessed}")
        if why:
            ev += f";low_because=depth_factor_fits;support={why}"
        elif short:
            ev += ";low_because=short"
        if copies:
            ev += f";other_copies={copies}"
            if r.copy_unique:
                ev += f";unique_at_copies={r.copy_unique};copies_llr={r.copy_llr}"
        rows.append(
            {
                "contig": r.contig, "start": int(r.start), "end": int(r.end), "class": cls, "role": role, "confidence": conf,
                "score": float(r.llr), "channel": "S1", "summary": what, "evidence": ev,
            }
        )
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


def segments_at_joins(segs: pd.DataFrame, joins: pd.DataFrame, pad: int = 5000) -> set:
    """Index of the CN_OR_DEPTH segments within `pad` of a join that no long read links."""
    keep = set()
    if not len(joins) or "verdict" not in joins.columns:
        return keep
    j = joins[joins["verdict"] == "unsupported"]
    by = {c: g[["gap_start", "gap_end"]].values for c, g in j.groupby("contig")}
    for i, r in segs[segs["kind"] == "CN_OR_DEPTH"].iterrows():
        e = by.get(r["contig"])
        if e is not None and np.any((e[:, 0] < r["end"] + pad) & (e[:, 1] > r["start"] - pad)):
            keep.add(i)
    return keep


def corroborate(segs: pd.DataFrame, others: list[pd.DataFrame], pad: int = 1000) -> set:
    """Index of the CN_OR_DEPTH segments that overlap an issue of another channel (a missing
    allele, haplotype or junction; a phase switch or misplaced piece; a false join)."""
    keep = set()
    parts = [t for t in others if t is not None and len(t)]
    if not parts:
        return keep
    ev = pd.concat(parts, ignore_index=True)
    by = {c: g[["start", "end"]].values for c, g in ev.groupby("contig")}
    for i, r in segs[segs["kind"] == "CN_OR_DEPTH"].iterrows():
        e = by.get(r["contig"])
        if e is not None and np.any((e[:, 0] < r["end"] + pad) & (e[:, 1] > r["start"] - pad)):
            keep.add(i)
    return keep


def depth_anomalies(segs: pd.DataFrame) -> pd.DataFrame:
    """Segments whose depth departs from the release without fitting an integer copy change."""
    keep = segs[segs["kind"].isin(["DEPTH", "IRREGULAR"])]  # ABSENT segments are issues, not anomalies
    cols = ["contig", "start", "end", "direction", "kind", "factor", "llr", "dominant_m", "copies_dominant", "classes", "n_kmers"]
    return keep[cols].reset_index(drop=True)


def issues_from_hapmers(switches: pd.DataFrame, pieces: pd.DataFrame) -> pd.DataFrame:
    rows = []
    whole = set()
    if len(pieces):
        for r in pieces[pieces["call"] == "other_parent"].itertuples(index=False):
            whole.add((r.contig, int(r.start), int(r.end)))
            rows.append(
                {
                    "contig": r.contig, "start": int(r.start), "end": int(r.end), "class": "MISPLACED_PIECE", "role": "locus",
                    "confidence": "high", "score": float(r.sites_other_parent), "channel": "S3",
                    "summary": f"scaffold piece carries the other parent's hapmers ({r.sites_other_parent} sites against "
                               f"{r.sites_scaffold_parent}) in a {r.scaffold_parent} scaffold",
                    "evidence": f"piece={r.piece};sites_other={r.sites_other_parent};sites_scaffold={r.sites_scaffold_parent}",
                }
            )
    if len(switches):
        for r in switches.itertuples(index=False):
            # a switch that is just a misplaced piece is reported once, as the piece
            if any(c == r.contig and s <= r.start and r.end <= e for c, s, e in whole):
                continue
            rows.append(
                {
                    "contig": r.contig, "start": int(r.start), "end": int(r.end), "class": "PHASE_SWITCH", "role": "locus",
                    "confidence": "high" if r.sites_segment_parent >= 10 else "low", "score": float(r.llr), "channel": "S3",
                    "summary": f"{r.sites_segment_parent} hapmer sites of {r.segment_parent} against {r.sites_scaffold_parent} "
                               f"of {r.scaffold_parent} in a {r.scaffold_parent} sequence",
                    "evidence": f"sites_segment_parent={r.sites_segment_parent};sites_scaffold_parent={r.sites_scaffold_parent}",
                }
            )
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


def issues_from_missing(clusters: pd.DataFrame, sample: int, run_gap: int = 10_000, run_min: int = 3,
                        prefer: pd.DataFrame | None = None, small: int = 50) -> pd.DataFrame:
    """Issues from the missing-k-mer channel.

    Runs of at least `run_min` small missing-sequence events less than `run_gap` apart are
    reported once, as a stretch of haplotype the release does not hold. An event whose flanks
    exist on both haplotypes is anchored on both; one placement is the `locus` (the one that
    overlaps an issue in `prefer`, else the first by name) and the other a `twin`.
    """
    rows = []
    if clusters is None or not len(clusters):
        return pd.DataFrame(columns=ISSUE_COLUMNS[1:])
    placed = clusters[clusters["contig"].notna()].copy()
    in_run = set()
    local = placed[(placed["kind"] == "MISSING_ALLELE") | ((placed["kind"] == "MISSING_SEQUENCE") & (placed["tracked_kmers"] < small))]
    local = local.sort_values(["contig", "start"])
    for contig, g in local.groupby("contig", sort=False):
        starts = g["start"].values
        brk = np.flatnonzero(np.diff(starts) > run_gap) + 1
        for idx in np.split(np.arange(len(g)), brk):
            if len(idx) < run_min:
                continue
            sub = g.iloc[idx]
            in_run.update(sub["cluster"].tolist())
            ev = f"alleles={len(sub)};reads={int(sub['n_reads_locus'].sum())};tracked_kmers={int(sub['tracked_kmers'].sum())};" \
                 f"median_count={float(sub['median_count'].median()):g}"
            if "frac_in_parent1" in sub and sub["frac_in_parent1"].notna().any():
                ev += f";in_parent1={sub['frac_in_parent1'].mean():.2f};in_parent2={sub['frac_in_parent2'].mean():.2f}"
            rows.append(
                {
                    "contig": contig, "start": int(sub["start"].min()), "end": int(sub["end"].max()), "class": "MISSING_HAPLOTYPE",
                    "role": "locus", "confidence": "high" if len(sub) >= 5 else "low", "score": float(len(sub)), "channel": "S2",
                    "summary": f"{len(sub)} alleles carried by the reads are absent from the release over "
                               f"{(int(sub['end'].max()) - int(sub['start'].min())) / 1000:.1f} kb: one haplotype of this stretch is not represented",
                    "evidence": ev, "_twins": frozenset(sub["twin_group"].tolist()),
                }
            )
    for r in placed.itertuples(index=False):
        if r.cluster in in_run:
            continue
        ev = f"reads={r.n_reads_locus};tracked_kmers={r.tracked_kmers};median_count={r.median_count:g}"
        if not pd.isna(getattr(r, "frac_in_parent1", np.nan)):
            ev += f";in_parent1={r.frac_in_parent1};in_parent2={r.frac_in_parent2}"
        if r.kind == "MISSING_JUNCTION":
            what = f"reads join this position to {r.contig2}:{int(r.pos2)}; no sequence in the release holds the junction"
            ev += f";mate={r.contig2}:{int(r.pos2)}"
        elif r.kind == "MISSING_ALLELE":
            size = int(r.size)
            what = "reads carry an allele here that no sequence in the release holds"
            if size:
                what += f" (the reads have {abs(size)} bases {'fewer' if size > 0 else 'more'} between the flanks)"
        else:
            what = f"reads run from here into about {int(r.tracked_kmers) * sample} k-mers of sequence found nowhere in the release"
        if isinstance(r.other_loci, str) and r.other_loci:
            ev += f";also_anchors={r.other_loci}"
        strong = r.n_reads_locus >= 6 and r.tracked_kmers >= 5
        conf = "high" if strong else "low"
        if r.kind == "MISSING_SEQUENCE":
            # Reads anchored on one side only are the weakest evidence of this channel: with few
            # reads nine in ten are artefacts (and all of them on a truth genome), so those stay
            # in the cluster table. The others are firm only when a parent's reads hold the
            # sequence. One exception among the weak ones needs the parents: recurrent artefacts
            # recur in the parents' reads too, so k-mers that are solid in the child and in neither
            # parent are few (some 45 events in a genome) and are where a de novo junction shows.
            p1, p2 = float(getattr(r, "frac_in_parent1", np.nan)), float(getattr(r, "frac_in_parent2", np.nan))
            in_neither = p1 == 0.0 and p2 == 0.0
            if not strong and not in_neither:
                continue
            conf = "high" if strong and max(p1, p2) >= 0.5 else "low"
            if in_neither:
                what += "; these k-mers are in neither parent's reads (a de novo junction or a change in the cell line)"
        rows.append(
            {
                "contig": r.contig, "start": int(max(r.start, 0)), "end": int(max(r.end, r.start + 1)), "class": r.kind,
                "role": "locus", "confidence": conf,
                "score": float(r.n_reads_locus), "channel": "S2", "summary": what, "evidence": ev,
                "_twins": frozenset([r.twin_group]),
            }
        )
    # one placement per event is the locus; the same event anchored on the other haplotype is its twin
    pref = {}
    if prefer is not None and len(prefer):
        pref = {c: g[["start", "end"]].values for c, g in prefer.groupby("contig")}

    def preferred(row) -> bool:
        e = pref.get(row["contig"])
        return e is not None and bool(np.any((e[:, 0] < row["end"] + 1000) & (e[:, 1] > row["start"] - 1000)))

    order = sorted(range(len(rows)), key=lambda i: (not preferred(rows[i]), rows[i]["contig"], rows[i]["start"]))
    taken: dict = {}
    for i in order:
        row = rows[i]
        owners = {taken[t] for t in row["_twins"] if t in taken}
        if owners and len([t for t in row["_twins"] if t in taken]) >= 0.5 * len(row["_twins"]):
            row["role"] = "twin"
            j = next(iter(owners))
            row["evidence"] += f";twin_of={rows[j]['contig']}:{rows[j]['start']}"
        else:
            for t in row["_twins"]:
                taken.setdefault(t, i)
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


def issues_from_unsupported(unsupported: pd.DataFrame, clusters: pd.DataFrame | None, k: int, w: Windows, ok: np.ndarray,
                            min_flank: float = 0.85, near: int = 60, covered: pd.DataFrame | None = None,
                            false_joins: bool = False) -> pd.DataFrame:
    """Issues from runs of assembly k-mers that no read supports, inside normally covered sequence.

    One wrong base and one false join both leave up to k unsupported k-mers. They differ in what
    the reads show: at a wrong base the reads carry another allele at the same place (a missing
    allele from S2); at a false join the reads carry nothing new, because the true neighbours of
    both sides exist elsewhere in the release.

    Only the first kind is reported by default. On a truth genome about a thousand runs of
    haplotype-specific k-mers have no reads although their flanks are normally covered (short
    stretches that the sequencing chemistry does not read), so an unsupported run without a
    second line of evidence is not an issue; every run stays in `unsupported_clusters.tsv`.
    Runs inside an issue of `covered` (copy-number segments already reported) are not repeated.
    """
    rows = []
    if unsupported is None or not len(unsupported):
        return pd.DataFrame(columns=ISSUE_COLUMNS[1:])
    name_to_ci = {n: i for i, n in enumerate(w.contigs["name"])}
    allele_pos = {}
    if clusters is not None and len(clusters):
        al = clusters[clusters["kind"] == "MISSING_ALLELE"]
        allele_pos = {c: np.sort(g["start"].values) for c, g in al.groupby("contig")}
    done = {}
    if covered is not None and len(covered):
        done = {c: g[["start", "end"]].values for c, g in covered.groupby("contig")}
    for r in unsupported.itertuples(index=False):
        if r.unsupported < r.tracked or not (r.flank_copies >= min_flank) or r.flank_kmers < 10:
            continue  # interleaved with supported k-mers, or the surroundings are poorly covered too
        d = done.get(r.contig)
        if d is not None:
            ov = np.minimum(d[:, 1], r.end) - np.maximum(d[:, 0], r.start)
            if np.any(ov >= 0.8 * (r.end - r.start)):
                continue
        ci = name_to_ci[r.contig]
        if not ok[int(w.contigs.at[ci, "win_start"]) + r.start // w.W]:
            continue
        a = allele_pos.get(r.contig)
        has_allele = False
        if a is not None:
            i = np.searchsorted(a, r.start - near)
            has_allele = i < len(a) and a[i] <= r.end + near
        if has_allele:
            cls, what = "BASE_ERROR", "assembly k-mers without read support, and the reads carry another allele here"
        elif false_joins:
            cls, what = "FALSE_JOIN", "no read supports the assembly's sequence across this point, and the reads carry no other allele here"
        else:
            continue
        rows.append(
            {
                "contig": r.contig, "start": int(r.start), "end": int(r.end), "class": cls, "role": "locus",
                "confidence": "high" if r.unsupported >= 6 else "low", "score": float(r.unsupported), "channel": "S1",
                "summary": what,
                "evidence": f"unsupported_tracked_kmers={r.unsupported};flank_copies={r.flank_copies};flank_kmers={r.flank_kmers}",
            }
        )
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[1:])


def combine(*tables: pd.DataFrame) -> pd.DataFrame:
    parts = [t for t in tables if t is not None and len(t)]
    if not parts:
        return pd.DataFrame(columns=ISSUE_COLUMNS)
    df = pd.concat(parts, ignore_index=True)
    df = df.sort_values(["contig", "start", "end", "class"]).reset_index(drop=True)
    df.insert(0, "id", [f"cp{i + 1:06d}" for i in range(len(df))])
    return df[ISSUE_COLUMNS]


def write_bed(issues: pd.DataFrame, path) -> None:
    """BED9, coloured by class. Columns 1-4 (contig, start, end, label) are what scorers read."""
    with open(path, "w") as f:
        f.write('track name=copyproof description="copyproof issues" itemRgb=On\n')
        for contig, start, end, cls, score in zip(issues["contig"], issues["start"], issues["end"], issues["class"], issues["score"]):
            sc = int(max(0, min(1000, round(float(score) * 10))))
            f.write(f"{contig}\t{start}\t{end}\t{cls}\t{sc}\t.\t{start}\t{end}\t{COLOURS.get(cls, '0,0,0')}\n")


def mask_runs(w: Windows, ok: np.ndarray) -> pd.DataFrame:
    """Runs of windows the copy-number channel could not assess."""
    rows = []
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        bad = ~ok[sl]
        if not bad.any():
            continue
        L = int(w.contigs.at[ci, "length"])
        cut = np.flatnonzero(np.diff(bad.astype(np.int8))) + 1
        starts = np.r_[0, cut]
        ends = np.r_[cut, len(bad)]
        name = w.contigs.at[ci, "name"]
        for s, e in zip(starts, ends):
            if bad[s]:
                rows.append((name, int(s) * w.W, min(int(e) * w.W, L)))
    return pd.DataFrame(rows, columns=["contig", "start", "end"])


def kmer_support(spectrum: pd.DataFrame, k: int, sample: int) -> dict:
    """Assembly k-mers without read support, and the consensus quality they imply (as Merqury does:
    a k-mer of the assembly that the reads lack marks at least one wrong base in k).
    Counts of 1 are stored as 0 by the engine, so at low depth this is a little pessimistic."""
    total = float(spectrum["n_kmers"].sum())
    unsupported = float(spectrum.loc[spectrum["count"] == 0, "n_kmers"].sum())
    if total == 0:
        return {}
    p_ok = (1.0 - unsupported / total) ** (1.0 / k)
    qv = float("inf") if p_ok >= 1 else -10.0 * np.log10(1.0 - p_ok)
    return {
        "distinct_kmers_estimated": int(total) * sample,
        "unsupported_kmers_estimated": int(unsupported) * sample,
        "unsupported_fraction": unsupported / total,
        "qv": round(qv, 2),
    }
