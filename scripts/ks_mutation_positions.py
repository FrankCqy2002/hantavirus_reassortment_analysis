#!/usr/bin/env python3
"""Mutation-position KS tests vs Uniform along each segment (all pairs + rank).

Combines mutation_positions_KS_all_pairs.py and mutation_positions_KS_rerank_readme.py.
H0: substitution positions among comparable sites are Uniform(0, n_comparable).
"""

from __future__ import annotations

import argparse
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

PKG_ROOT = Path(__file__).resolve().parents[1]


def resolve_data_dir(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    env = os.environ.get("ANDV_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (PKG_ROOT / ".." / "hantavirus" / "Alignment" / "analyses" / "06_v5_label_aln").resolve()


def short_name(r: pd.Series) -> str:
    d = str(r["sample_descriptor"])
    if "Chile-9717869" in str(r["full_label"]):
        return "Chile-9717869"
    if str(r.get("sequence_group", "")).startswith("Cruise"):
        return "cruise"
    if d.startswith("ANDV "):
        return d.replace("ANDV ", "")
    return d


def mut_comp_coords(seq_a: str, seq_b: str):
    a, b = seq_a.upper(), seq_b.upper()
    mut = []
    n_comp = 0
    for x, y in zip(a, b):
        if x in "-N" or y in "-N":
            continue
        if x != y:
            mut.append(n_comp)
        n_comp += 1
    return np.asarray(mut, dtype=np.float64), n_comp


def ks_uniform(mut: np.ndarray, n_comp: int):
    n = len(mut)
    if n < 2 or n_comp <= 0:
        return n, np.nan, np.nan, (n / n_comp if n_comp else np.nan)
    D, p = stats.kstest(mut, "uniform", args=(0.0, float(n_comp)))
    return n, float(D), float(p), n / n_comp


def load_scope(piet: Path, scope: str):
    if scope == "full":
        attr = pd.read_csv(piet / "metadata" / "andv_attributes_v5.csv")
        alns = {
            seg: {
                r.id: str(r.seq).upper()
                for r in SeqIO.parse(piet / "alignments" / f"{seg}_v5_label_aln.fasta", "fasta")
            }
            for seg in "SML"
        }
    else:
        base = piet / "clade_III"
        attr = pd.read_csv(base / "metadata" / "clade_III_attributes_v5.csv")
        alns = {
            seg: {
                r.id: str(r.seq).upper()
                for r in SeqIO.parse(
                    base / "alignments" / f"{seg}_cladeIII_v5.fasta", "fasta"
                )
            }
            for seg in "SML"
        }
    return attr, alns


def compute_all_pairs(attr: pd.DataFrame, alns: dict[str, dict[str, str]]) -> pd.DataFrame:
    tip_info = {}
    for _, r in attr.iterrows():
        lab = short_name(r)
        tip_info[lab] = {
            "S": r["meta_S"],
            "M": r["meta_M"],
            "L": r["meta_L"],
            "clade": r.get("andv_clade", ""),
            "full_label": r["full_label"],
        }

    labels = sorted(tip_info.keys())
    print(f"n tips: {len(labels)}")
    seqs = {seg: {lab: alns[seg][tip_info[lab][seg]] for lab in labels} for seg in "SML"}

    rows = []
    n_pairs = len(labels) * (len(labels) - 1) // 2
    print(f"Computing {n_pairs} pairs × 3 segments...")
    for i, (a, b) in enumerate(combinations(labels, 2)):
        row = {
            "seq_a": a,
            "seq_b": b,
            "clade_a": tip_info[a]["clade"],
            "clade_b": tip_info[b]["clade"],
        }
        for seg in ("S", "M", "L"):
            mut, n_comp = mut_comp_coords(seqs[seg][a], seqs[seg][b])
            n_mut, D, p, pdist = ks_uniform(mut, n_comp)
            row[f"{seg}_n_mut"] = n_mut
            row[f"{seg}_n_comparable"] = n_comp
            row[f"{seg}_p_distance"] = pdist
            row[f"{seg}_ks_D"] = D
            row[f"{seg}_ks_p"] = p
        rows.append(row)
        if (i + 1) % 500 == 0:
            print(f"  {i + 1}/{n_pairs}")
    return pd.DataFrame(rows)


def rerank(out: pd.DataFrame, min_mut: int) -> pd.DataFrame:
    out = out.copy()
    out["L_neglog10_p"] = -np.log10(out["L_ks_p"].clip(lower=1e-300))
    out["S_neglog10_p"] = -np.log10(out["S_ks_p"].clip(lower=1e-300))
    out["M_neglog10_p"] = -np.log10(out["M_ks_p"].clip(lower=1e-300))

    eligible = out["L_n_mut"] >= min_mut
    ranked = (
        out.loc[eligible]
        .sort_values(["L_ks_D", "L_neglog10_p"], ascending=[False, False])
        .copy()
    )
    ranked.insert(0, "rank_L", np.arange(1, len(ranked) + 1))
    rest = out.loc[~eligible].copy()
    rest.insert(0, "rank_L", np.nan)
    final = pd.concat([ranked, rest], ignore_index=True)

    cols = [
        "rank_L",
        "seq_a",
        "seq_b",
        "clade_a",
        "clade_b",
        "S_n_mut",
        "S_n_comparable",
        "S_p_distance",
        "S_ks_D",
        "S_ks_p",
        "S_neglog10_p",
        "M_n_mut",
        "M_n_comparable",
        "M_p_distance",
        "M_ks_D",
        "M_ks_p",
        "M_neglog10_p",
        "L_n_mut",
        "L_n_comparable",
        "L_p_distance",
        "L_ks_D",
        "L_ks_p",
        "L_neglog10_p",
    ]
    cols = [c for c in cols if c in final.columns]
    return final[cols], int(eligible.sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--min-mut", type=int, default=20, help="Min L mutations for rank_L")
    ap.add_argument("--scope", choices=["full", "cladeIII"], default="full")
    args = ap.parse_args()

    piet = resolve_data_dir(args.data_dir)
    outdir = (args.out_dir or (PKG_ROOT / "results" / "tables")).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    attr, alns = load_scope(piet, args.scope)
    raw = compute_all_pairs(attr, alns)
    final, n_elig = rerank(raw, args.min_mut)

    out_path = outdir / "mutation_positions_KS_vs_uniform_all_pairs.csv"
    final.to_csv(out_path, index=False, float_format="%.6g")
    print(f"\nWrote {out_path} ({len(final)} pairs)")

    readme = outdir / "mutation_positions_KS_README.txt"
    readme.write_text(
        "Mutation-position KS test vs Uniform along comparable sites.\n"
        "H0: differences are uniformly distributed along the segment.\n"
        "Larger ks_D / smaller ks_p => more clustered (non-uniform) mutations.\n"
        f"rank_L: among pairs with L_n_mut >= {args.min_mut}, sorted by L_ks_D descending.\n"
        "Pairs with fewer L mutations have blank rank_L (KS unreliable at tiny n).\n"
    )
    print(f"Wrote {readme}")
    print(f"Eligible for rank_L (L_n_mut>={args.min_mut}): {n_elig} / {len(final)}")
    print("\nTop 15 by rank_L:")
    top = final.dropna(subset=["rank_L"]).head(15)
    cols = ["rank_L", "seq_a", "seq_b", "clade_a", "clade_b", "L_n_mut", "L_ks_D", "L_ks_p", "L_p_distance"]
    print(top[[c for c in cols if c in top.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
