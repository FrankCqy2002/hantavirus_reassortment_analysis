#!/usr/bin/env python3
"""Sliding-window p-distance profiles for any isolate pair (S/M/L).

Generalized from clean_f615_line425_allpair_sw.py and plot_p1236_cruise_sliding_window.py.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

PKG_ROOT = Path(__file__).resolve().parents[1]
COLORS = {"S": "#2B6CB0", "M": "#6B7C3A", "L": "#C0392B"}


def resolve_data_dir(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    env = os.environ.get("ANDV_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (PKG_ROOT / ".." / "hantavirus" / "Alignment" / "analyses" / "06_v5_label_aln").resolve()


def resolve_sample(attr: pd.DataFrame, key: str) -> pd.Series:
    key = key.strip()
    hit = attr.loc[attr["sample_descriptor"].astype(str) == key]
    if len(hit) == 1:
        return hit.iloc[0]
    if key.lower() in {"cruise", "cruise ship", "cruise-ship", "outbreak"}:
        hit = attr.loc[
            attr["sequence_group"].astype(str).str.contains("Cruise", case=False, na=False)
        ]
        if len(hit) == 1:
            return hit.iloc[0]
    hit = attr.loc[
        attr["sample_descriptor"].astype(str).str.contains(re.escape(key), case=False, na=False)
        | attr["full_label"].astype(str).str.contains(re.escape(key), case=False, na=False)
        | attr["meta_S"].astype(str).str.contains(re.escape(key), case=False, na=False)
        | attr["meta_M"].astype(str).str.contains(re.escape(key), case=False, na=False)
        | attr["meta_L"].astype(str).str.contains(re.escape(key), case=False, na=False)
    ]
    if len(hit) == 1:
        return hit.iloc[0]
    if len(hit) > 1:
        sub = hit.loc[hit["full_label"].astype(str).str.contains(re.escape(key), case=False, na=False)]
        if len(sub) == 1:
            return sub.iloc[0]
    raise SystemExit(
        f"Could not uniquely resolve {key!r}: {len(hit)} hits\n"
        f"{hit[['sample_descriptor', 'full_label']].to_string() if len(hit) else ''}"
    )


def short_label(r: pd.Series) -> str:
    d = str(r["sample_descriptor"])
    if d and d != "nan":
        if d.startswith("ANDV "):
            return d.replace("ANDV ", "")
        return d
    fl = str(r["full_label"])
    m = re.search(r"Chile[-_]?9717869", fl, re.I)
    if m:
        return "Chile-9717869"
    if str(r.get("sequence_group", "")).startswith("Cruise"):
        return "cruise"
    return fl


def window_pdists(seq_a: str, seq_b: str, win: int, step: int) -> pd.DataFrame:
    a, b = seq_a.upper(), seq_b.upper()
    comparable = [i for i in range(len(a)) if a[i] not in "-N" and b[i] not in "-N"]
    rows = []
    for start in range(0, len(comparable) - win + 1, step):
        idxs = comparable[start : start + win]
        diffs = sum(1 for i in idxs if a[i] != b[i])
        rows.append(
            {
                "aln_start": idxs[0] + 1,
                "aln_end": idxs[-1] + 1,
                "aln_mid": (idxs[0] + idxs[-1]) / 2 + 1,
                "mutations": diffs,
                "p_distance": diffs / win,
            }
        )
    return pd.DataFrame(rows)


def slug(label: str) -> str:
    return re.sub(r"[^\w.\-]+", "-", label).strip("-")


def overall_pdist(seq_a: str, seq_b: str) -> float:
    a, b = seq_a.upper(), seq_b.upper()
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


def plot_publication(fig_path: Path, lab_a: str, lab_b: str, panels: list[tuple[str, pd.DataFrame, float]]) -> None:
    """Three stacked panels matching results/sliding_window publication PNGs."""
    plt.rcParams["font.family"] = "Arial"
    fig, axes = plt.subplots(3, 1, figsize=(3.53, 2.61), sharex=False)
    ymax = max(float(df.p_distance.max()) for _, df, _ in panels)
    ylim = max(0.08, ymax * 1.15)
    for ax, (seg, df, overall) in zip(axes, panels):
        color = COLORS[seg]
        ax.fill_between(df.aln_mid, df.p_distance, color=color, alpha=0.12, linewidth=0)
        ax.plot(df.aln_mid, df.p_distance, color=color, lw=1.0, solid_capstyle="butt")
        ax.axhline(overall, color="#888888", ls="--", lw=0.7)
        ax.set_ylim(0, ylim)
        ax.set_yticks([0.0, 0.05] if ylim <= 0.10 else [0.0, 0.05, 0.10])
        ax.tick_params(axis="both", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(
            0.02,
            0.90,
            f"{seg} segment (dist = {overall:.4f})",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
    axes[0].set_title(f"{lab_a} vs {lab_b}", fontsize=9, pad=2)
    axes[1].set_ylabel("p-distance", fontsize=10)
    axes[2].set_xlabel("Alignment position (window midpoint)", fontsize=10)
    fig.tight_layout(pad=0.35)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    if fig_path.suffix.lower() == ".png":
        fig.savefig(fig_path.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seq-a", required=True, help="Sample descriptor / alias for sequence A")
    ap.add_argument("--seq-b", required=True, help="Sample descriptor / alias for sequence B")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--scope", choices=["full", "cladeIII"], default="full")
    ap.add_argument(
        "--per-segment",
        action="store_true",
        help="Also write single-segment PNG/PDF panels",
    )
    ap.add_argument(
        "--publication",
        action="store_true",
        help="Write the stacked publication panel used in results/sliding_window",
    )
    args = ap.parse_args()

    piet = resolve_data_dir(args.data_dir)
    outdir = Path(args.out_dir) if args.out_dir else (PKG_ROOT / "results" / "sliding_window")
    outdir.mkdir(parents=True, exist_ok=True)

    if args.scope == "full":
        attr = pd.read_csv(piet / "metadata" / "andv_attributes_v5.csv")
        alns = {
            seg: {r.id: str(r.seq) for r in SeqIO.parse(piet / "alignments" / f"{seg}_v5_label_aln.fasta", "fasta")}
            for seg in "SML"
        }
    else:
        base = piet / "clade_III"
        attr = pd.read_csv(base / "metadata" / "clade_III_attributes_v5.csv")
        alns = {
            seg: {
                r.id: str(r.seq)
                for r in SeqIO.parse(
                    base / "alignments" / f"{seg}_cladeIII_v5.fasta", "fasta"
                )
            }
            for seg in "SML"
        }

    ra = resolve_sample(attr, args.seq_a)
    rb = resolve_sample(attr, args.seq_b)
    lab_a, lab_b = short_label(ra), short_label(rb)
    print(f"A: {lab_a} clade={ra.get('andv_clade', '')} tips={{S:{ra['meta_S']}, M:{ra['meta_M']}, L:{ra['meta_L']}}}")
    print(f"B: {lab_b} clade={rb.get('andv_clade', '')} tips={{S:{rb['meta_S']}, M:{rb['meta_M']}, L:{rb['meta_L']}}}")

    all_rows = []
    overalls: dict[str, float] = {}
    fig, axes = plt.subplots(3, 1, figsize=(12, 9))
    for ax, seg in zip(axes, ("S", "M", "L")):
        tip_a, tip_b = ra[f"meta_{seg}"], rb[f"meta_{seg}"]
        if tip_a not in alns[seg] or tip_b not in alns[seg]:
            raise SystemExit(f"Missing tip(s) in {seg} alignment: {tip_a!r}, {tip_b!r}")
        df = window_pdists(alns[seg][tip_a], alns[seg][tip_b], args.window, args.step)
        overalls[seg] = overall_pdist(alns[seg][tip_a], alns[seg][tip_b])
        df.insert(0, "segment", seg)
        if args.publication:
            df.insert(1, "seq_a", lab_a)
            df.insert(2, "seq_b", lab_b)
            df.insert(3, "overall_p_distance", overalls[seg])
        all_rows.append(df)
        mean, mx = df.p_distance.mean(), df.p_distance.max()
        print(f"  {seg}: mean={mean:.4f}, max={mx:.4f}, windows={len(df)}")

        ax.plot(df.aln_mid, df.p_distance, color=COLORS[seg], lw=1.3)
        ax.fill_between(df.aln_mid, df.p_distance, color=COLORS[seg], alpha=0.15)
        ax.axhline(mean, color="#888888", ls="--", lw=0.8)
        ax.set_ylabel("p-distance")
        ax.set_xlabel(f"Alignment position (midpoint of {args.window}-bp window)")
        ax.set_title(f"Segment {seg}: {lab_a} vs {lab_b}")
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if args.per_segment:
            fig_s, ax_s = plt.subplots(figsize=(12, 3.2))
            ax_s.plot(df.aln_mid, df.p_distance, color=COLORS[seg], lw=1.3)
            ax_s.fill_between(df.aln_mid, df.p_distance, color=COLORS[seg], alpha=0.15)
            ax_s.set_xlabel(f"Alignment position (midpoint of {args.window}-bp window)")
            ax_s.set_ylabel("p-distance")
            ax_s.set_title(f"Segment {seg}: {lab_a} vs {lab_b}")
            ax_s.set_ylim(bottom=0)
            ax_s.spines["top"].set_visible(False)
            ax_s.spines["right"].set_visible(False)
            fig_s.tight_layout()
            stem_s = f"{slug(lab_a)}_vs_{slug(lab_b)}_sliding_window_pdist_{seg}"
            fig_s.savefig(outdir / f"{stem_s}.png", dpi=300, bbox_inches="tight")
            plt.close(fig_s)

    stem = f"{slug(lab_a)}_vs_{slug(lab_b)}_sliding_window_pdist"
    if args.publication:
        stem = f"{stem}_publication"
        plt.close(fig)
        plot_publication(
            outdir / f"{stem}.png",
            lab_a,
            lab_b,
            [(seg, df, overalls[seg]) for seg, df in zip(("S", "M", "L"), all_rows)],
        )
    else:
        clade_a = ra.get("andv_clade", "")
        clade_b = rb.get("andv_clade", "")
        fig.suptitle(
            f"Sliding-window p-distance: {lab_a} (Clade {clade_a}) vs {lab_b} (Clade {clade_b})  "
            f"({args.window} bp, step {args.step})",
            fontsize=12,
            y=1.01,
        )
        fig.tight_layout()
        fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
    pd.concat(all_rows, ignore_index=True).to_csv(outdir / f"{stem}.csv", index=False)
    print(f"Wrote {outdir / stem}.png / .csv")


if __name__ == "__main__":
    main()
