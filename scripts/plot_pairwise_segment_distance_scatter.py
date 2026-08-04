#!/usr/bin/env python3
"""All-pairs S/M/L segment p-distance scatters with binomial λ testing.

For every tip pair, compute observed p-distance / mutation counts on S, M, and L
(ignoring gap/N sites), fit a conditional binomial rate-ratio λ, and flag
BH-significant outliers for publication-style scatters (y = λx; red points =
BH q < alpha).
"""

from __future__ import annotations

import argparse
import os
import re
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib.lines import Line2D
from scipy.optimize import minimize_scalar
from scipy.stats import norm

PKG_ROOT = Path(__file__).resolve().parents[1]
COMPARISONS = [
    ("S", "M"),
    ("S", "L"),
    ("M", "L"),
]
AXIS_LABELS = {
    "S": "S-segment pairwise distance",
    "M": "M-segment pairwise distance",
    "L": "L-segment pairwise distance",
}


def resolve_data_dir(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    env = os.environ.get("ANDV_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (
        PKG_ROOT / ".." / "hantavirus" / "Alignment" / "analyses" / "06_v5_label_aln"
    ).resolve()


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
                for r in SeqIO.parse(base / "alignments" / f"{seg}_cladeIII_v5.fasta", "fasta")
            }
            for seg in ("S", "M", "L")
        }
    return meta, alignments


def resolve_keys(meta: pd.DataFrame, keys: list[str]) -> set[str]:
    """Resolve substring keys to sample_descriptor values."""
    hit_desc: set[str] = set()
    for key in keys:
        key = key.strip()
        if not key:
            continue
        hit = meta.loc[
            meta["sample_descriptor"].astype(str).str.contains(re.escape(key), case=False, na=False)
            | meta["full_label"].astype(str).str.contains(re.escape(key), case=False, na=False)
        ]
        hit_desc.update(hit["sample_descriptor"].astype(str).tolist())
    return hit_desc


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


def annotate_involves(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["involves_p1236"] = out["descriptor_a"].astype(str).str.contains("1236", case=False) | out[
        "descriptor_b"
    ].astype(str).str.contains("1236", case=False)
    cruise = out["group_a"].astype(str).str.contains("Cruise", case=False) | out[
        "group_b"
    ].astype(str).str.contains("Cruise", case=False)
    out["involves_cruise"] = cruise & ~out["involves_p1236"]
    return out


def bh_fdr(pvals: np.ndarray, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=bool), np.array([])
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out <= alpha, out


def ols_lambda(px, py) -> float:
    x = np.asarray(px, float)
    y = np.asarray(py, float)
    denom = np.sum(x * x)
    if denom <= 0:
        return 1.0
    return float(np.sum(x * y) / denom)


def pi_of(lam: float, e_x: np.ndarray, e_y: np.ndarray) -> np.ndarray:
    return (lam * e_y) / (e_x + lam * e_y)


def mu_tpi(lam: float, y_x, y_y, e_x, e_y) -> np.ndarray:
    t = np.asarray(y_x, float) + np.asarray(y_y, float)
    return t * pi_of(lam, np.asarray(e_x, float), np.asarray(e_y, float))


def binom_mle_lambda(y_x, y_y, e_x, e_y) -> tuple[float, float]:
    """Conditional binomial MLE for common λ across pairs."""
    y_x = np.asarray(y_x, float)
    y_y = np.asarray(y_y, float)
    e_x = np.asarray(e_x, float)
    e_y = np.asarray(e_y, float)

    def nll(lam: float) -> float:
        if lam <= 0:
            return 1e300
        pi = np.clip(pi_of(lam, e_x, e_y), 1e-12, 1 - 1e-12)
        return -float(np.sum(y_y * np.log(pi) + y_x * np.log(1 - pi)))

    res = minimize_scalar(nll, bounds=(1e-4, 50.0), method="bounded", options={"xatol": 1e-10})
    return float(res.x), float(-res.fun)


def binomial_test_panel(
    df: pd.DataFrame,
    xseg: str,
    yseg: str,
    mu_min: float,
    fdr_alpha: float,
) -> pd.DataFrame:
    """Fit λ on all pairs with Tπ > mu_min; test each pair; return enriched table."""
    mx, my = f"mut_{xseg}", f"mut_{yseg}"
    nx, ny = f"n_{xseg}", f"n_{yseg}"
    px, py = f"p_{xseg}", f"p_{yseg}"

    lam0 = ols_lambda(df[px], df[py])
    mu0 = mu_tpi(lam0, df[mx], df[my], df[nx], df[ny])
    fit0 = df.loc[mu0 > mu_min]
    if fit0.empty:
        raise SystemExit(f"No pairs with Tπ>{mu_min:g} for preliminary λ fit ({xseg} vs {yseg})")

    lam, _ = binom_mle_lambda(fit0[mx], fit0[my], fit0[nx], fit0[ny])
    mu1 = mu_tpi(lam, df[mx], df[my], df[nx], df[ny])
    fit = df.loc[mu1 > mu_min]
    lam, _ = binom_mle_lambda(fit[mx], fit[my], fit[nx], fit[ny])
    lam_ols = ols_lambda(fit[px], fit[py])
    n_fit = int(len(fit))

    rows = []
    for _, r in df.iterrows():
        y_x = float(r[mx])
        y_y = float(r[my])
        e_x = float(r[nx])
        e_y = float(r[ny])
        t = y_x + y_y
        pi = pi_of(lam, np.array([e_x]), np.array([e_y]))[0] if (e_x + lam * e_y) > 0 else np.nan
        mu = t * pi if np.isfinite(pi) else np.nan
        var = t * pi * (1 - pi) if np.isfinite(pi) else np.nan
        if np.isfinite(var) and var > 0:
            z = (y_y - t * pi) / np.sqrt(var)
            pval = float(2 * norm.cdf(-abs(z)))
        else:
            z = np.nan
            pval = 1.0
        used = np.isfinite(mu) and mu > mu_min
        rows.append(
            {
                "seq_a": r["seq_a"],
                "seq_b": r["seq_b"],
                "descriptor_a": r["descriptor_a"],
                "descriptor_b": r["descriptor_b"],
                "hhpc_a": r["hhpc_a"],
                "hhpc_b": r["hhpc_b"],
                px: r[px],
                py: r[py],
                mx: y_x,
                my: y_y,
                nx: e_x,
                ny: e_y,
                "T": t,
                "pi": pi,
                "involves_p1236": bool(r["involves_p1236"]),
                "involves_cruise": bool(r["involves_cruise"]),
                "lambda": lam,
                "lambda_ols_ref": lam_ols,
                "mu_Y": mu,
                "testable_mu_gt_threshold": bool(used),
                "testable_mu_gt_20": bool(used),  # legacy column name
                "used_in_lambda_fit": bool(used),
                "z_binomial": z,
                "pvalue_raw": pval,
            }
        )
    out = pd.DataFrame(rows)

    out["pvalue"] = np.nan
    out["pvalue_BH"] = np.nan
    out["sig_0.05"] = False
    mask = out["testable_mu_gt_threshold"].to_numpy()
    if mask.any():
        reject, q = bh_fdr(out.loc[mask, "pvalue_raw"].to_numpy(), alpha=fdr_alpha)
        out.loc[mask, "pvalue"] = out.loc[mask, "pvalue_raw"].to_numpy()
        out.loc[mask, "pvalue_BH"] = q
        out.loc[mask, "sig_0.05"] = reject

    out = out.sort_values(
        ["sig_0.05", "pvalue"], ascending=[False, True], na_position="last"
    ).reset_index(drop=True)

    n_testable = int(mask.sum())
    n_sig = int(out["sig_0.05"].sum())
    print(
        f"{xseg} vs {yseg}: λ̂={lam:.6f} (OLS ref={lam_ols:.6f}); "
        f"n_fit={n_fit}; sig={n_sig}/{n_testable} at BH q<{fdr_alpha:g} (Tπ>{mu_min:g})"
    )
    return out


def plot_binomial_publication(
    out: pd.DataFrame,
    xseg: str,
    yseg: str,
    outfile: Path,
    mu_min: float,
    fdr_alpha: float,
) -> None:
    """Publication-style panel: grey points, red BH-significant, y=λx."""
    plot = out[out["mu_Y"] > mu_min].copy()
    if plot.empty:
        print(f"Skip plot {xseg} vs {yseg}: no pairs with Tπ>{mu_min:g}")
        return

    lam = float(plot["lambda"].iloc[0])
    color_other = "#B0B0B0"
    color_sig = "#C0392B"
    color_line = "#1a1a1a"

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    ns = plot[~plot["sig_0.05"]]
    sig = plot[plot["sig_0.05"]]
    ax.scatter(
        ns[f"p_{xseg}"],
        ns[f"p_{yseg}"],
        s=14,
        c=color_other,
        alpha=0.85,
        linewidths=0,
        zorder=2,
    )
    if len(sig):
        ax.scatter(
            sig[f"p_{xseg}"],
            sig[f"p_{yseg}"],
            s=22,
            c=color_sig,
            edgecolors="black",
            linewidths=0.7,
            zorder=5,
        )

    xmax = float(np.nanmax(plot[f"p_{xseg}"]))
    ymax = float(np.nanmax(plot[f"p_{yseg}"]))
    lim_x = xmax * 1.06 if xmax > 0 else 0.01
    lim_y = ymax * 1.08 if ymax > 0 else 0.01
    xs = np.linspace(0, lim_x, 200)
    ax.plot(xs, lam * xs, color=color_line, lw=1.1, ls="--", zorder=3)

    ax.set_xlabel(AXIS_LABELS[xseg])
    ax.set_ylabel(AXIS_LABELS[yseg])
    ax.set_xlim(0, lim_x)
    ax.set_ylim(0, lim_y)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(
        0.98,
        0.04,
        f"λ = {lam:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=color_line,
    )
    if plot["sig_0.05"].any():
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color_sig,
                markeredgecolor="black",
                markersize=7,
                markeredgewidth=0.9,
                label=f"BH q < {fdr_alpha:g}",
            )
        ]
        ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=8)

    fig.tight_layout(pad=0.3)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outfile} (n={len(plot)}, n_sig={int(plot['sig_0.05'].sum())})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="06_v5_label_aln directory (default: ANDV_DATA_DIR or ../hantavirus/Alignment/...)",
    )
    ap.add_argument(
        "--scope",
        choices=["full", "cladeIII"],
        default="cladeIII",
        help="Tip set (default: cladeIII)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: results/segment_distance_scatter)",
    )
    ap.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Sample substrings to drop before pairwise comparisons",
    )
    ap.add_argument(
        "--mu-min",
        type=float,
        default=20.0,
        help="Minimum expected Y count Tπ for λ fit and BH testing (default: 20)",
    )
    ap.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="Benjamini–Hochberg FDR threshold for significance (default: 0.05)",
    )
    args = ap.parse_args()

    piet = resolve_data_dir(args.data_dir)
    out_dir = (args.out_dir or (PKG_ROOT / "results" / "segment_distance_scatter")).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    meta, alignments = load_scope(piet, args.scope)
    if args.exclude:
        drop = resolve_keys(meta, args.exclude)
        if drop:
            meta = meta[~meta["sample_descriptor"].astype(str).isin(drop)].copy()
            print(f"Excluded {len(drop)} tip(s): {sorted(drop)}")

    pairs = annotate_involves(compute_pairs(meta, alignments))
    stem_prefix = "CURATED" if args.scope == "full" else "CladeIII"
    csv_path = out_dir / f"{stem_prefix}_pairwise_segment_pdist.csv"
    pairs.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} (n_pairs={len(pairs)}, n_tips={len(meta)})")

    print(f"Binomial test: mu_min={args.mu_min:g}, fdr_alpha={args.fdr_alpha:g}")
    for xseg, yseg in COMPARISONS:
        out = binomial_test_panel(
            pairs,
            xseg,
            yseg,
            mu_min=args.mu_min,
            fdr_alpha=args.fdr_alpha,
        )
        test_csv = out_dir / f"{stem_prefix}_{xseg}_vs_{yseg}_lambda_binomial_test.csv"
        out.to_csv(test_csv, index=False)
        print(f"Wrote {test_csv}")
        plot_binomial_publication(
            out,
            xseg,
            yseg,
            out_dir / f"{stem_prefix}_pairwise_{xseg}_vs_{yseg}_binomial_publication.png",
            mu_min=args.mu_min,
            fdr_alpha=args.fdr_alpha,
        )


if __name__ == "__main__":
    main()
