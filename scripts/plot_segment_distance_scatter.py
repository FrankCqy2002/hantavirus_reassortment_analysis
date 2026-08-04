#!/usr/bin/env python3
"""Scatter plots of pairwise distances to a reference across S/M/L segments.

Supports the full curated database (default) or the Clade III subset.
Data root: --data-dir, env ANDV_DATA_DIR, or ../Alingments Piet relative to package.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib.lines import Line2D

PKG_ROOT = Path(__file__).resolve().parents[1]

GROUP_COLORS = {
    "Chilean sequences": "#6B7C3A",
    "Argentinean sequences": "#1F4E79",
    "Cruise ship outbreak sequences": "#C0392B",
}

GROUP_SHORT = {
    "Chilean sequences": "Chilean",
    "Argentinean sequences": "Argentinean",
    "Cruise ship outbreak sequences": "Cruise ship outbreak",
}

# NCBI RefSeq assembly GCF_000850405.1 = Chile-9717869 (NC_003466/7/8)
NCBI_REF_ALIASES = {
    "gcf_000850405.1",
    "gcf_000850405",
    "gcf000850405.1",
    "nc_003466",
    "nc_003466.1",
    "nc_003467",
    "nc_003467.2",
    "nc_003468",
    "nc_003468.2",
}


def resolve_data_dir(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    env = os.environ.get("ANDV_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (PKG_ROOT / ".." / "Alingments Piet").resolve()


def short_label(row: pd.Series) -> str:
    desc = str(row["sample_descriptor"])
    if desc.startswith("ANDV "):
        return desc.replace("ANDV ", "")
    if desc.startswith("Switzerland"):
        return "Switzerland-Case 7"
    return desc


def p_distance(seq_a: str, seq_b: str) -> float:
    """Observed substitutions/site, ignoring gap/N sites in either sequence."""
    a = seq_a.upper()
    b = seq_b.upper()
    if len(a) != len(b):
        raise ValueError("Aligned sequences must have equal length")
    valid = diffs = 0
    for x, y in zip(a, b):
        if x in "-N" or y in "-N":
            continue
        valid += 1
        if x != y:
            diffs += 1
    if valid == 0:
        return float("nan")
    return diffs / valid


def load_alignment(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def parse_mldist(path: Path) -> pd.DataFrame:
    lines = path.read_text().strip().splitlines()
    n = int(lines[0].strip())
    names: list[str] = []
    mat = np.zeros((n, n), dtype=float)
    for i, line in enumerate(lines[1 : n + 1]):
        parts = line.split()
        names.append(parts[0])
        mat[i, :] = [float(x) for x in parts[1 : n + 1]]
    return pd.DataFrame(mat, index=names, columns=names)


def identical_tip_map(alignment: dict[str, str], present: set[str]) -> dict[str, str]:
    """Map tips missing from mldist to an identical sequence that is present."""
    by_seq: dict[str, list[str]] = {}
    for tip, seq in alignment.items():
        by_seq.setdefault(seq, []).append(tip)
    mapping: dict[str, str] = {}
    for tips in by_seq.values():
        reps = [t for t in tips if t in present]
        if not reps:
            continue
        rep = reps[0]
        for t in tips:
            mapping[t] = rep
    return mapping


def resolve_reference(meta: pd.DataFrame, ref_key: str) -> pd.Series:
    key = ref_key.strip()
    if key.lower().replace(" ", "") in NCBI_REF_ALIASES or key.upper().startswith("GCF_000850405"):
        key = "Chile-9717869"
    hit = meta.loc[meta["sample_descriptor"].astype(str) == key]
    if len(hit) == 1:
        return hit.iloc[0]
    if key.lower() in {"cruise", "cruise ship", "cruise-ship", "outbreak", "switzerland"}:
        hit = meta.loc[
            meta["sequence_group"].astype(str).str.contains("Cruise ship", case=False, na=False)
        ]
        if len(hit) == 1:
            return hit.iloc[0]
    hit = meta.loc[
        meta["sample_descriptor"].astype(str).str.contains(key, case=False, na=False)
        | meta["full_label"].astype(str).str.contains(key, case=False, na=False)
    ]
    if len(hit) == 1:
        return hit.iloc[0]
    raise SystemExit(
        f"Expected exactly one reference for {ref_key!r}, got {len(hit)}: "
        f"{hit['sample_descriptor'].tolist() if len(hit) else []}"
    )


def mldist_lookup(
    mldist: pd.DataFrame,
    tip_map: dict[str, str],
    ref_tip: str,
    tip: str,
) -> float:
    a = tip_map.get(ref_tip, ref_tip)
    b = tip_map.get(tip, tip)
    if a not in mldist.index or b not in mldist.columns:
        return float("nan")
    return float(mldist.loc[a, b])


def distances_to_reference(
    meta: pd.DataFrame,
    alignments: dict[str, dict[str, str]],
    mldists: dict[str, pd.DataFrame],
    tip_maps: dict[str, dict[str, str]],
    ref_key: str,
) -> tuple[pd.DataFrame, pd.Series]:
    ref = resolve_reference(meta, ref_key)

    rows = []
    for _, row in meta.iterrows():
        if row["sample_descriptor"] == ref["sample_descriptor"]:
            continue
        out = {
            "sample_descriptor": row["sample_descriptor"],
            "label": short_label(row),
            "full_label": row["full_label"],
            "sequence_group": row["sequence_group"],
            "highlighted_HHPC_cluster": bool(row["highlighted_HHPC_cluster"]),
            "year_group": row["year_group"],
            "psi_code": row["psi_code"],
            "andv_clade": row.get("andv_clade", ""),
            "reference": short_label(ref),
            "reference_ncbi": "GCF_000850405.1",
        }
        for seg in ("S", "M", "L"):
            tip = row[f"meta_{seg}"]
            ref_tip = ref[f"meta_{seg}"]
            out[f"pdist_{seg}"] = p_distance(alignments[seg][tip], alignments[seg][ref_tip])
            out[f"mldist_{seg}"] = mldist_lookup(mldists[seg], tip_maps[seg], ref_tip, tip)
        rows.append(out)
    return pd.DataFrame(rows), ref


def add_panel(
    ax,
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    xlabel: str,
    ylabel: str,
    annotate: str,
    zoom_to_data: bool = False,
):
    xmin, xmax = float(df[xcol].min()), float(df[xcol].max())
    ymin, ymax = float(df[ycol].min()), float(df[ycol].max())
    if zoom_to_data:
        lo = min(xmin, ymin)
        hi = max(xmax, ymax)
        pad = max((hi - lo) * 0.18, 0.001)
        lim_lo, lim_hi = lo - pad, hi + pad
    else:
        lim_lo, lim_hi = 0.0, max(xmax, ymax) * 1.12

    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color="#888888", lw=1, ls="--", zorder=0)
    for _, r in df.iterrows():
        color = GROUP_COLORS.get(r["sequence_group"], "#333333")
        marker = "D" if r["highlighted_HHPC_cluster"] else "o"
        edge = "black" if r["highlighted_HHPC_cluster"] else color
        ax.scatter(
            r[xcol],
            r[ycol],
            c=color,
            marker=marker,
            s=55 if r["highlighted_HHPC_cluster"] else 45,
            edgecolors=edge,
            linewidths=0.8,
            zorder=2,
        )
        do_ann = annotate == "all" or (annotate == "hhpc" and r["highlighted_HHPC_cluster"])
        if do_ann:
            ax.annotate(
                r["label"],
                (r[xcol], r[ycol]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=5.5 if annotate == "all" else 6.5,
                color="#222222",
            )
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_scatters(
    df: pd.DataFrame,
    out: Path,
    metric: str,
    title: str,
    annotate: str = "hhpc",
    zoom_to_data: bool = False,
):
    prefix = "pdist" if metric == "p-distance" else "mldist"
    unit = "p-distance (subst./site)" if metric == "p-distance" else "ML distance (subst./site)"
    pairs = [
        (f"{prefix}_S", f"{prefix}_M", f"S {unit}", f"M {unit}", "S vs M"),
        (f"{prefix}_M", f"{prefix}_L", f"M {unit}", f"L {unit}", "M vs L"),
        (f"{prefix}_S", f"{prefix}_L", f"S {unit}", f"L {unit}", "S vs L"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.6))
    for ax, (xc, yc, xl, yl, panel) in zip(axes, pairs):
        add_panel(ax, df, xc, yc, xl, yl, annotate=annotate, zoom_to_data=zoom_to_data)
        ax.set_title(panel, fontsize=11)

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8, label=GROUP_SHORT[g])
        for g, c in GROUP_COLORS.items()
    ]
    legend_handles += [
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
        Line2D([0], [0], color="#888888", lw=1, ls="--", label="y = x"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
        fontsize=8.5,
    )
    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")
    print(f"Wrote {out.with_suffix('.pdf')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Alingments Piet directory (default: ANDV_DATA_DIR or ../Alingments Piet)",
    )
    ap.add_argument(
        "--scope",
        choices=["full", "cladeIII"],
        default="full",
        help="full curated database (default) or Clade III subset",
    )
    ap.add_argument("--metric", choices=["p-distance", "ml"], default="p-distance")
    ap.add_argument(
        "--reference",
        default="GCF_000850405.1",
        help="Reference isolate (GCF_000850405.1 / Chile-9717869, sample_descriptor, or 'cruise')",
    )
    ap.add_argument(
        "--hhpc-only",
        action="store_true",
        help="Restrict scatter points to highlighted HHPC cluster isolates",
    )
    ap.add_argument(
        "--annotate",
        choices=["none", "hhpc", "all"],
        default=None,
        help="Point labels (default: hhpc for full, all for cladeIII)",
    )
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    piet = resolve_data_dir(args.data_dir)

    if args.scope == "full":
        meta = pd.read_csv(piet / "andv_attributes_joined_metadata.csv")
        alignments = {
            seg: load_alignment(piet / f"{seg}_complete_aln_CURATED.fasta") for seg in ("S", "M", "L")
        }
        mldists = {
            seg: parse_mldist(piet / "phylogeny" / "iqtree" / seg / f"{seg}_ML.mldist")
            for seg in ("S", "M", "L")
        }
        out_dir = args.out_dir or (PKG_ROOT / "results" / "segment_distance")
        scope_label = "ANDV curated"
        stem_prefix = "CURATED"
    else:
        base = piet / "clade_III"
        meta = pd.read_csv(base / "metadata" / "clade_III_attributes_joined.csv")
        alignments = {
            seg: load_alignment(base / "alignments" / f"{seg}_cladeIII_subset_from_CURATED.fasta")
            for seg in ("S", "M", "L")
        }
        mldists = {
            seg: parse_mldist(base / "iqtree" / seg / f"{seg}_ML.mldist") for seg in ("S", "M", "L")
        }
        out_dir = args.out_dir or (PKG_ROOT / "results" / "segment_distance")
        scope_label = "ANDV Clade III"
        stem_prefix = "CladeIII"

    tip_maps = {
        seg: identical_tip_map(alignments[seg], set(mldists[seg].index)) for seg in ("S", "M", "L")
    }

    df, ref = distances_to_reference(
        meta, alignments, mldists, tip_maps, ref_key=args.reference
    )
    ref_label = short_label(ref)
    if str(args.reference).upper().startswith("GCF_000850405") or "9717869" in ref_label:
        ref_slug = "Chile-9717869_GCF_000850405.1"
    elif "Cruise" in str(ref["sequence_group"]):
        ref_slug = "cruise"
    else:
        ref_slug = ref_label.replace("/", "-").replace(" ", "_")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"{ref_slug}_segment_pairwise_distances.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} (reference={ref_label}; n={len(df)} vs {len(meta)} curated)")

    if args.hhpc_only:
        df = df[df["highlighted_HHPC_cluster"]].copy()
        hhpc_csv = out_dir / f"{ref_slug}_HHPC_segment_pairwise_distances.csv"
        df.to_csv(hhpc_csv, index=False)
        print(f"Wrote {hhpc_csv}")

    annotate = args.annotate
    if annotate is None:
        annotate = "all" if args.scope == "cladeIII" or args.hhpc_only else "hhpc"

    metric = args.metric
    dist_word = "pairwise distance" if metric == "p-distance" else "ML distance"
    if args.hhpc_only:
        title = f"{scope_label} HHPC cluster — {dist_word} to {ref_label} across segments"
        stem = f"{stem_prefix}_HHPC_{ref_slug}_segment_distance_scatter"
    else:
        title = f"{scope_label} — {dist_word} to {ref_label} (GCF_000850405.1) across segments"
        if "9717869" not in ref_label and not str(args.reference).upper().startswith("GCF"):
            title = f"{scope_label} — {dist_word} to {ref_label} across segments"
        stem = f"{stem_prefix}_{ref_slug}_segment_distance_scatter"
    if metric == "ml":
        stem += "_ML"
    plot_scatters(
        df,
        out_dir / f"{stem}.png",
        metric,
        title,
        annotate=annotate,
        zoom_to_data=args.hhpc_only,
    )


if __name__ == "__main__":
    main()
