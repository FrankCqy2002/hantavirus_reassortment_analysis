# Hantavirus (ANDV) reassortment analysis

Self-contained publication scripts for Andes virus (ANDV) segment-distance
scatters, sliding-window p-distances, and Kolmogorov–Smirnov (KS) tests that
probe non-uniform divergence along S/M/L segments (reassortment / mosaic signal).

## Environment

```bash
conda env create -f environment.yml
conda activate hantavirus_reassortment
# or: pip install -r requirements.txt
```

## Data layout

Point `--data-dir` at the **Alingments Piet** directory (or set `ANDV_DATA_DIR`).
Default: `../Alingments Piet` next to this package. See `data/README.md`.

Outputs go under `--out-dir` (default: `results/` in this package).

## Methods (brief)

### Segment distance scatter
For each non-reference isolate, compute pairwise distances to a reference on
S, M, and L (observed p-distance ignoring gap/N sites; optional IQ-TREE ML
distances from `.mldist`). Plot S–M, M–L, and S–L scatters with a y=x line.
Points far from the diagonal indicate segment-incongruent distances.

### Sliding-window p-distance
Along comparable (non-gap/N) sites, slide a window (default 100 nt, step 10),
count differences / window size, and plot the profile for S/M/L for any pair
resolved from metadata (`--seq-a`, `--seq-b`).

### Sliding-window KS test
For each `*_sliding_window_pdist.csv`, test whether window p-distances follow
**Uniform(0, max)** (primary) or Uniform(min, max) (secondary). Larger KS *D*
implies stronger departure from a flat window-distance profile.

### Mutation-position KS test
For all tip pairs, collect 0-based coordinates of substitutions among
comparable sites and test those positions against **Uniform(0, n_comparable)**.
Larger *D* / smaller *p* suggests clustered mutations. Pairs are ranked by L
segment KS *D*, requiring at least `--min-mut` L mutations (default 20).

## Usage

```bash
export ANDV_DATA_DIR="/path/to/Alingments Piet"   # optional

# Distance scatters (full curated DB or Clade III)
python scripts/plot_segment_distance_scatter.py --scope full --metric p-distance \
  --data-dir "$ANDV_DATA_DIR" --out-dir results/segment_distance

python scripts/plot_segment_distance_scatter.py --scope cladeIII --metric ml \
  --reference Chile-9717869 --data-dir "$ANDV_DATA_DIR" --out-dir results/segment_distance

# Sliding window for any pair
python scripts/plot_sliding_window_pdist.py --seq-a Chile-9717869 --seq-b p1236 \
  --window 100 --step 10 --scope full --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/sliding_window

# KS on sliding-window CSVs
python scripts/ks_sliding_window.py \
  --input-dir results/sliding_window --out-dir results/tables

# Mutation-position KS (all pairs + rank)
python scripts/ks_mutation_positions.py --scope full --min-mut 20 \
  --data-dir "$ANDV_DATA_DIR" --out-dir results/tables
```

## Repository contents

| Path | Role |
|------|------|
| `scripts/plot_segment_distance_scatter.py` | S/M/L distance scatters to a reference |
| `scripts/plot_sliding_window_pdist.py` | Pairwise sliding-window p-distance plots |
| `scripts/ks_sliding_window.py` | KS of window p-distances vs Uniform |
| `scripts/ks_mutation_positions.py` | KS of mutation positions vs Uniform |
| `data/` | Placeholder; real inputs via `--data-dir` |
| `results/` | Default output (gitignored) |
