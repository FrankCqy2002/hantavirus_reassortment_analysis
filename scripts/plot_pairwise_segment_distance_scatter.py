#!/usr/bin/env python3
"""All-pairs S/M/L segment p-distance scatters (not distance-to-reference).

For every tip pair in the chosen scope, compute observed p-distance on S, M,
and L (ignoring gap/N sites), then plot S–M, M–L, and S–L scatters. Points far
from y=x indicate segment-incongruent distances (possible reassortment signal).
"""

from __future__ import annotations

import argparse
import os
import re
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from matplotlib.lines import Line2D

PKG_ROOT = Path(__file__).resolve().parents[1]


def resolve_data_dir(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    env = os.environ.get("ANDV_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (PKG_ROOT / ".." / "hantavirus" / "Alignment" / "analyses" / "06_v5_label_aln").resolve()


def short_label(row: pd.Series) -> str:
    desc = str(row["sample_descriptor"])
    if desc.startswith("ANDV "):
        return desc.replace("ANDV ", "")
    if desc.startswith("Switzerland"):
        return "Switzerland-Case 7"
    if "Chile-9717869" in str(row.get("full_label", "")):
        return "Chile-9717869"
    if str(row.get("sequence_group", "")).startswith("Cruise"):
        return "cruise"
    return desc


def p_distance(seq_a: str, seq_b: str) -> tuple[float, int, int]:
    """Return (p-distance, n_mutations, n_comparable), ignoring gap/N sites."""
    a, b = seq_a.upper(), seq_b.upper()
    valid = diffs = 0
    for x, y in zip(a, b):
        if x in "-N" or y in "-N":
            continue
        valid += 1
        if x != y:
            diffs += 1
    if valid == 0:
        return float("nan"), 0, 0
    return diffs / valid, diffs, valid


def load_scope(piet: Path, scope: str) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    if scope == "full":
        meta = pd.read_csv(piet / "metadata" / "andv_attributes_v5.csv")
        alignments = {
            seg: {
                r.id: str(r.seq).upper()
                for r in SeqIO.parse(piet / "alignments" / f"{seg}_v5_label_aln.fasta", "fasta")
            }
            for seg in ("S", "M", "L")
        }
    else:
        base = piet / "clade_III"
        meta = pd.read_csv(base / "metadata" / "clade_III_attributes_v5.csv")
        alignments = {
            seg: {
                r.id: str(r.seq).upper()
                for r in SeqIO.parse(
                    base / "alignments" / f"{seg}_cladeIII_v5.fasta", "fasta"
                )
            }
            for seg in ("S", "M", "L")
        }
    return meta, alignments


def resolve_exclude(meta: pd.DataFrame, keys: list[str]) -> set[str]:
    """Return sample_descriptor values to drop."""
    drop: set[str] = set()
    for key in keys:
        key = key.strip()
        if not key:
            continue
        hit = meta.loc[
            meta["sample_descriptor"].astype(str).str.contains(re.escape(key), case=False, na=False)
            | meta["full_label"].astype(str).str.contains(re.escape(key), case=False, na=False)
        ]
        drop.update(hit["sample_descriptor"].astype(str).tolist())
    return drop


def compute_pairs(meta: pd.DataFrame, alignments: dict[str, dict[str, str]]) -> pd.DataFrame:
    rows = []
    for (_, a), (_, b) in combinations(meta.iterrows(), 2):
        row = {
            "seq_a": short_label(a),
            "seq_b": short_label(b),
            "descriptor_a": a["sample_descriptor"],
            "descriptor_b": b["sample_descriptor"],
            "clade_a": a.get("andv_clade", ""),
            "clade_b": b.get("andv_clade", ""),
            "group_a": a.get("sequence_group", ""),
            "group_b": b.get("sequence_group", ""),
            "hhpc_a": bool(a.get("highlighted_HHPC_cluster", False)),
            "hhpc_b": bool(b.get("highlighted_HHPC_cluster", False)),
            "same_clade": a.get("andv_clade", "") == b.get("andv_clade", ""),
        }
        for seg in ("S", "M", "L"):
            tip_a, tip_b = a[f"meta_{seg}"], b[f"meta_{seg}"]
            p, n_mut, n_comp = p_distance(alignments[seg][tip_a], alignments[seg][tip_b])
            row[f"p_{seg}"] = p
            row[f"mut_{seg}"] = n_mut
            row[f"n_{seg}"] = n_comp
        rows.append(row)
    return pd.DataFrame(rows)


def pair_style(r: pd.Series, highlight: list[str]) -> tuple[str, str, bool]:
    """Return (color, marker, edge_black) for a pair row."""
    involves_cruise = ("Cruise" in str(r["group_a"])) or ("Cruise" in str(r["group_b"]))
    involves_highlight = False
    for key in highlight:
        key_l = key.lower()
        if key_l in str(r["seq_a"]).lower() or key_l in str(r["seq_b"]).lower():
            involves_highlight = True
            break
        if key_l in str(r["descriptor_a"]).lower() or key_l in str(r["descriptor_b"]).lower():
            involves_highlight = True
            break

    if involves_cruise:
        color = "#C0392B"
    elif involves_highlight:
        color = "#2B6CB0"
    else:
        color = "#6B7C3A"
    marker = "D" if (r["hhpc_a"] or r["hhpc_b"]) else "o"
    return color, marker, involves_highlight


def plot_panel(ax, pairs: pd.DataFrame, xseg: str, yseg: str, highlight: list[str], title: str):
    xmax = max(pairs[f"p_{xseg}"].max(), pairs[f"p_{yseg}"].max()) * 1.12
    if not (xmax > 0):
        xmax = 0.01
    ax.plot([0, xmax], [0, xmax], color="#888888", lw=1, ls="--", zorder=0)
    for _, r in pairs.iterrows():
        color, marker, edge_black = pair_style(r, highlight)
        ax.scatter(
            r[f"p_{xseg}"],
            r[f"p_{yseg}"],
            c=color,
            marker=marker,
            s=40,
            edgecolors="black" if edge_black else color,
            linewidths=0.6,
            alpha=0.85,
            zorder=2,
        )
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, xmax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{xseg} p-distance")
    ax.set_ylabel(f"{yseg} p-distance")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def legend_handles(highlight: list[str]) -> list[Line2D]:
    hl = ", ".join(highlight) if highlight else "highlight"
    return [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#6B7C3A", markersize=8, label="Other pairs"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#2B6CB0",
            markeredgecolor="black",
            markersize=8,
            label=f"Involves {hl}",
        ),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#C0392B", markersize=8, label="Involves cruise"),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor="#888888",
            markeredgecolor="black",
            markersize=8,
            label="HHPC cluster",
        ),
        Line2D([0], [0], color="#888888", ls="--", label="y = x"),
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="06_v5_label_aln directory (default: ANDV_DATA_DIR or ../hantavirus/Alignment/analyses/06_v5_label_aln)",
    )
    ap.add_argument(
        "--scope",
        choices=["full", "cladeIII"],
        default="cladeIII",
        help="Tip set: full curated DB or Clade III subset (default: cladeIII)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: results/segment_distance_scatter)",
    )
    ap.add_argument(
        "--highlight",
        nargs="*",
        default=["p1236"],
        help="Sample substrings to highlight in blue (default: p1236). Cruise is always red.",
    )
    ap.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Sample substrings to exclude from pairwise comparisons",
    )
    ap.add_argument(
        "--combined-only",
        action="store_true",
        help="Write only the 3-panel combined figure (skip per-panel PNGs)",
    )
    args = ap.parse_args()

    piet = resolve_data_dir(args.data_dir)
    out_dir = (args.out_dir or (PKG_ROOT / "results" / "segment_distance_scatter")).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    meta, alignments = load_scope(piet, args.scope)
    if args.exclude:
        drop = resolve_exclude(meta, args.exclude)
        if drop:
            meta = meta[~meta["sample_descriptor"].astype(str).isin(drop)].copy()
            print(f"Excluded {len(drop)} tip(s): {sorted(drop)}")

    pairs = compute_pairs(meta, alignments)
    stem_prefix = "CURATED" if args.scope == "full" else "CladeIII"
    csv_path = out_dir / f"{stem_prefix}_pairwise_segment_pdist.csv"
    pairs.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} (n_pairs={len(pairs)}, n_tips={len(meta)})")

    scope_label = "ANDV curated" if args.scope == "full" else "ANDV Clade III"
    highlight = list(args.highlight or [])

    if not args.combined_only:
        for xseg, yseg in [("S", "M"), ("M", "L"), ("S", "L")]:
            fig, ax = plt.subplots(figsize=(5.8, 5.4))
            plot_panel(ax, pairs, xseg, yseg, highlight, f"{scope_label} pairwise {xseg} vs {yseg}")
            ax.legend(handles=legend_handles(highlight), frameon=False, fontsize=8, loc="upper left")
            fig.tight_layout()
            stem = f"{stem_prefix}_pairwise_{xseg}_vs_{yseg}_scatter"
            fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"Wrote {out_dir / stem}.png")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))
    for ax, (xseg, yseg) in zip(axes, [("S", "M"), ("M", "L"), ("S", "L")]):
        plot_panel(ax, pairs, xseg, yseg, highlight, f"{xseg} vs {yseg}")
    fig.legend(
        handles=legend_handles(highlight),
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.05),
        fontsize=8.5,
    )
    fig.suptitle(f"{scope_label} pairwise segment p-distances", y=1.02)
    fig.tight_layout()
    stem = f"{stem_prefix}_pairwise_S_M_L_scatter"
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_dir / stem}.png")


if __name__ == "__main__":
    main()
