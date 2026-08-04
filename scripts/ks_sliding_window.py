#!/usr/bin/env python3
"""KS tests of sliding-window p-distances vs Uniform(0, max).

Recovered from f615_line438_SW_KS_SW_KS2.py. Reads *_sliding_window_pdist.csv
files produced by plot_sliding_window_pdist.py.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PKG_ROOT = Path(__file__).resolve().parents[1]


def ks_vs_uniform(x: np.ndarray) -> dict:
    """KS vs Uniform(0, max(x)); also Uniform(min, max)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return dict(
            ks_D=np.nan,
            ks_p=np.nan,
            ks_D_minmax=np.nan,
            ks_p_minmax=np.nan,
            n=len(x),
            mean=np.nan,
            max=np.nan,
        )
    xmax = float(np.max(x))
    xmin = float(np.min(x))
    if xmax <= 0:
        D0, p0 = 0.0, 1.0
    else:
        D0, p0 = stats.kstest(x, "uniform", args=(0.0, xmax))
    span = xmax - xmin
    if span <= 0:
        D1, p1 = 0.0, 1.0
    else:
        D1, p1 = stats.kstest(x, "uniform", args=(xmin, span))
    return dict(
        ks_D=float(D0),
        ks_p=float(p0),
        ks_D_minmax=float(D1),
        ks_p_minmax=float(p1),
        n=int(len(x)),
        mean=float(np.mean(x)),
        max=xmax,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing *_sliding_window_pdist.csv files",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for KS tables (default: results/tables)",
    )
    args = ap.parse_args()

    sw = args.input_dir.expanduser().resolve()
    outdir = (args.out_dir or (PKG_ROOT / "results" / "tables")).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    csvs = sorted(sw.glob("*_sliding_window_pdist.csv"))
    if not csvs:
        raise SystemExit(f"No *_sliding_window_pdist.csv under {sw}")

    rows = []
    for csv in csvs:
        stem = csv.name.replace("_sliding_window_pdist.csv", "")
        m = re.match(r"(.+)_vs_(.+)$", stem)
        if not m:
            pair_a, pair_b = stem, ""
        else:
            pair_a, pair_b = m.group(1), m.group(2)
        df = pd.read_csv(csv)
        row = {"pair": stem, "seq_a": pair_a, "seq_b": pair_b}
        for seg in ("S", "M", "L"):
            sub = df.loc[df["segment"] == seg, "p_distance"].to_numpy()
            res = ks_vs_uniform(sub)
            row[f"{seg}_ks_D"] = res["ks_D"]
            row[f"{seg}_ks_p"] = res["ks_p"]
            row[f"{seg}_mean_p"] = res["mean"]
            row[f"{seg}_max_p"] = res["max"]
            row[f"{seg}_n_windows"] = res["n"]
            row[f"{seg}_ks_D_minmax"] = res["ks_D_minmax"]
            row[f"{seg}_ks_p_minmax"] = res["ks_p_minmax"]
        rows.append(row)

    out = pd.DataFrame(rows)
    out = out.sort_values("L_ks_D", ascending=False).reset_index(drop=True)
    out.insert(0, "rank_L", np.arange(1, len(out) + 1))

    main_tbl = out[
        [
            "rank_L",
            "pair",
            "seq_a",
            "seq_b",
            "S_ks_D",
            "S_ks_p",
            "M_ks_D",
            "M_ks_p",
            "L_ks_D",
            "L_ks_p",
            "S_mean_p",
            "M_mean_p",
            "L_mean_p",
        ]
    ].copy()

    out_path = outdir / "sliding_window_KS_vs_uniform.csv"
    full_path = outdir / "sliding_window_KS_vs_uniform_full.csv"
    main_tbl.to_csv(out_path, index=False, float_format="%.6g")
    out.to_csv(full_path, index=False, float_format="%.6g")
    print(f"Wrote {out_path}")
    print(f"Wrote {full_path}")
    print("\nTop 10 by L KS D (vs Uniform(0, max)):")
    print(main_tbl.head(10).to_string(index=False))
    print("\nNote: H0 = window p-distances ~ Uniform(0, observed max). Larger D => less uniform.")


if __name__ == "__main__":
    main()
