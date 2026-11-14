"""Evidence cards: one JSON record and one small figure per issue.

A card shows what a curator needs to judge an issue: how many copies the reads imply along the
locus for the k-mers the release holds once or twice, where each parent's hapmers fall, and where
reads carry sequence the release lacks. The figure is a set of small panels on one shared position
axis; the JSON next to it holds the same numbers.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .model import Model
from .windows import Windows

# Palette and chrome: the first three categorical slots (validated as a set) for the k-mer classes;
# everything else is drawn in neutral ink, so colour only ever means "k-mer class".
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
WASH = "#f0efec"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]


def class_colour(m: int) -> str:
    """Colour follows the k-mer class, not its rank in a card: two copies blue, one copy orange,
    anything else aqua."""
    return {2: SERIES[0], 1: SERIES[1]}.get(m, SERIES[2])


def class_style(m: int, drawn: list[int]) -> dict:
    """Line and marker of a class. The classes above two copies share aqua; when a card draws two
    of them, the higher one is dashed with open markers, so the two are told apart without a
    fourth hue."""
    higher = m > 2 and any(2 < other < m for other in drawn)
    if higher:
        return {"ls": (0, (4, 2)), "marker": "o", "mfc": SURFACE, "mec": class_colour(m), "mew": 1.0, "ms": 3.8}
    return {"ls": "-", "marker": "o", "mfc": class_colour(m), "mec": SURFACE, "mew": 0.7, "ms": 3.6}

EVENT_ROWS = [
    ("MISSING_JUNCTION", "missing junction"),
    ("MISSING_HAPLOTYPE", "missing haplotype"),
    ("MISSING_ALLELE", "missing allele"),
    ("MISSING_SEQUENCE", "missing sequence"),
    ("FALSE_JOIN", "unsupported join"),
    ("BASE_ERROR", "unsupported k-mers"),
    ("GAP", "gap (N)"),
]


def _view(start: int, end: int, length: int, W: int, max_view: int = 2_000_000) -> tuple[int, int]:
    span = max(end - start, W)
    pad = int(min(max(10_000, span), max_view // 2))
    a = max(0, (start - pad) // W * W)
    b = min(length, -(-(end + pad) // W) * W)
    return a, b


def card_data(issue: dict, w: Windows, model: Model, ok: np.ndarray, issues: pd.DataFrame, gaps: pd.DataFrame,
              max_points: int = 200) -> dict:
    """Numbers behind one card."""
    names = list(w.contigs["name"])
    ci = names.index(issue["contig"])
    L = int(w.contigs.at[ci, "length"])
    sl = w.slice_of(ci)
    a, b = _view(int(issue["start"]), int(issue["end"]), L, w.W)
    wa, wb = a // w.W, -(-b // w.W)
    step = max(1, -(-(wb - wa) // max_points))
    nb = -(-(wb - wa) // step)
    idx = np.arange(wa, wb)
    grp = (idx - wa) // step
    n = w.n[sl][wa:wb].astype(float)
    S = w.S[sl][wa:wb].astype(float)
    f = model.gc_factor(w.gc[sl][wa:wb].astype(float))
    okv = ok[sl][wa:wb]
    use = okv[:, None] & np.ones_like(n, dtype=bool)
    n_use = np.where(use, n, 0.0)
    # classes to draw: the commonest inside the issue, and the haplotype-unique class
    iw0, iw1 = int(issue["start"]) // w.W - wa, -(-int(issue["end"]) // w.W) - wa
    inside = n_use[max(iw0, 0):max(iw1, 1)].sum(axis=0)
    inside[0] = 0
    order = [int(j) for j in np.argsort(-inside) if inside[j] >= 20]
    classes = order[:1]
    if 1 in order and 1 not in classes:
        classes.append(1)
    elif len(order) > 1 and len(classes) < 2:
        classes.append(order[1])
    series = {}
    for j in classes:
        nj = np.bincount(grp, weights=n_use[:, j], minlength=nb)
        Sj = np.bincount(grp, weights=np.where(use[:, j], S[:, j], 0.0), minlength=nb)
        Ej = np.bincount(grp, weights=model.lam * f * n_use[:, j], minlength=nb)
        with np.errstate(invalid="ignore", divide="ignore"):
            series[j] = np.where(nj >= 8, Sj / Ej, np.nan)
    xs = a + (np.arange(nb) * step + step / 2.0) * w.W
    masked = np.bincount(grp, weights=(~okv).astype(float), minlength=nb) / np.bincount(grp, minlength=nb)
    hap = None
    if w.hap is not None:
        h = w.hap[sl][wa:wb].astype(float)
        hap = [np.bincount(grp, weights=h[:, 0], minlength=nb).tolist(), np.bincount(grp, weights=h[:, 1], minlength=nb).tolist()]
    ev = issues[(issues["contig"] == issue["contig"]) & (issues["end"] > a) & (issues["start"] < b)]
    events = {k: [] for k, _ in EVENT_ROWS}
    for cls, s0, e0 in zip(ev["class"], ev["start"], ev["end"]):
        if cls in events:
            events[cls].append([int(s0), int(e0)])
    g = gaps[(gaps["contig"] == issue["contig"]) & (gaps["end"] > a) & (gaps["start"] < b)] if len(gaps) else gaps
    events["GAP"] = [[int(s0), int(e0)] for s0, e0 in zip(g["start"], g["end"])] if len(g) else []
    return {
        "issue": {k: (v.item() if hasattr(v, "item") else v) for k, v in issue.items()},
        "view": [int(a), int(b)], "bin_bp": int(step * w.W), "x": xs.tolist(),
        "copies_by_multiplicity": {str(j): [None if np.isnan(v) else round(float(v), 3) for v in series[j]] for j in classes},
        "not_assessable_fraction": np.round(masked, 3).tolist(),
        "hapmer_sites": hap, "events": {k: v for k, v in events.items() if v},
        "depth_per_copy": round(model.lam, 3),
    }


def draw_card(card: dict, path, parent_names: tuple[str, str] | None = None) -> None:
    """Render a card as a PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"], "font.size": 9, "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    })
    iss = card["issue"]
    x = np.array(card["x"])
    a, b = card["view"]
    unit, div = ("Mb", 1e6) if b > 2_000_000 else ("kb", 1e3)
    has_hap = card["hapmer_sites"] is not None and parent_names is not None
    rows = [k for k, _ in EVENT_ROWS if k in card["events"]]
    heights = [3.0] + ([1.6] if has_hap else []) + ([0.28 * len(rows) + 0.35] if rows else [])
    fig, axes = plt.subplots(len(heights), 1, figsize=(9.2, 0.9 + sum(heights) * 0.72), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.16}, squeeze=False)
    axes = axes[:, 0]
    for ax in axes:
        ax.axvspan(iss["start"] / div, max(iss["end"], iss["start"] + (b - a) / 400) / div, color=WASH, lw=0, zorder=0)
        ax.tick_params(length=3, width=0.6)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.6)

    # 1. copies implied by the reads
    ax = axes[0]
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    top = 2.5
    drawn = [int(j) for j in card["copies_by_multiplicity"]]
    for i, (j, ys) in enumerate(card["copies_by_multiplicity"].items()):
        m = int(j)
        y = np.array([np.nan if v is None else v for v in ys], dtype=float)
        ax.axhline(m, color=AXIS, lw=0.9, zorder=1)
        ax.text(b / div, m, f" release: {m}", va="center", ha="left", fontsize=8, color=INK2, clip_on=False)
        ax.plot(x / div, y, color=class_colour(m), lw=1.3, zorder=3, label=f"k-mers the release holds {m}x",
                **class_style(m, drawn))
        if np.isfinite(y).any():
            top = max(top, float(np.nanmax(y)) * 1.12, m + 1.2)
    ax.set_ylim(0, min(top, 12))
    ax.set_ylabel("copies in the reads")
    na = np.array(card["not_assessable_fraction"])
    if (na > 0.5).any():
        w_bin = card["bin_bp"] / div
        ax.bar(x[na > 0.5] / div, 0.035 * ax.get_ylim()[1], width=w_bin, bottom=0, color=MUTED, lw=0, zorder=2)
        ax.text(a / div, 0.05 * ax.get_ylim()[1], " grey bar: not assessable", fontsize=7.5, color=MUTED, va="bottom", ha="left")
    if len(card["copies_by_multiplicity"]) > 1:
        ax.legend(frameon=False, fontsize=8, loc="upper left", labelcolor=INK2, handlelength=1.6, borderaxespad=0.2)
    elif card["copies_by_multiplicity"]:
        m = next(iter(card["copies_by_multiplicity"]))
        ax.text(0.005, 0.97, f"k-mers the release holds {m}x", transform=ax.transAxes, fontsize=8, color=INK2, va="top")
    k = 1
    # 2. parental hapmer sites, mirrored
    if has_hap:
        ax = axes[k]
        k += 1
        up, down = np.array(card["hapmer_sites"][0]), np.array(card["hapmer_sites"][1])
        wbar = card["bin_bp"] / div * 0.8
        ax.bar(x / div, up, width=wbar, color=INK2, lw=0)
        ax.bar(x / div, -down, width=wbar, color=INK2, lw=0)
        ax.axhline(0, color=AXIS, lw=0.9)
        lim = max(3.0, float(max(up.max(), down.max())) * 1.25)
        ax.set_ylim(-lim, lim)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, integer=True, symmetric=True))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}"))
        ax.yaxis.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_ylabel("hapmer sites")
        ax.text(0.005, 0.95, f"up: only in {parent_names[0]}", transform=ax.transAxes, fontsize=8, color=INK2, va="top")
        ax.text(0.005, 0.05, f"down: only in {parent_names[1]}", transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
    # 3. where reads and release disagree in sequence
    if rows:
        ax = axes[k]
        labels = dict(EVENT_ROWS)
        for r, key in enumerate(rows):
            for s0, e0 in card["events"][key]:
                if e0 - s0 > (b - a) / 150:
                    ax.plot([s0 / div, e0 / div], [r, r], color=INK2, lw=3.2, solid_capstyle="butt")
                else:
                    ax.plot([(s0 + e0) / 2 / div], [r], marker="|", ms=9, mew=1.3, color=INK2, ls="none")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([labels[r] for r in rows], fontsize=8, color=INK2)
        ax.set_ylim(len(rows) - 0.4, -0.6)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[-1].set_xlim(a / div, b / div)
    axes[-1].set_xlabel(f"{iss['contig']} position ({unit})")
    title = f"{iss.get('id', '')}  {iss['class']}  {iss['contig']}:{int(iss['start']):,}-{int(iss['end']):,}"
    fig.suptitle(title, x=0.01, ha="left", fontsize=10.5, fontweight="bold", y=0.995)
    summary = str(iss.get("summary", ""))
    if len(summary) > 150:
        summary = summary[:147] + "..."
    fig.text(0.01, 0.955 if len(heights) > 1 else 0.93, summary, fontsize=8.5, color=INK2, ha="left", va="top")
    fig.subplots_adjust(left=0.135 if rows else 0.085, right=0.9, top=0.86 if len(heights) > 1 else 0.8,
                        bottom=0.11 if len(heights) > 1 else 0.17)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_cards(issues: pd.DataFrame, w: Windows, model: Model, ok: np.ndarray, gaps: pd.DataFrame, outdir,
                parent_names: tuple[str, str] | None = None, limit: int | None = None, figures: bool = True,
                only: pd.DataFrame | None = None) -> int:
    """Write card JSON (and PNG) for the issues in `only` (default: all); every issue of the
    table is drawn as context. Returns how many cards were written."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for row in (issues if only is None else only).to_dict("records"):
        if limit is not None and n >= limit:
            break
        card = card_data(row, w, model, ok, issues, gaps)
        (out / f"{row['id']}.json").write_text(json.dumps(card))
        if figures:
            draw_card(card, out / f"{row['id']}.png", parent_names)
        n += 1
    return n
