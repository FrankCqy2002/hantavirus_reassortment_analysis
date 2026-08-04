# Hantavirus (ANDV) reassortment analysis

Python package to test potential reassortment in Andes virus (ANDV) using
pairwise segment p-distance scatters, sliding-window profiles, and KS tests.

## Environment

Requires [uv](https://docs.astral.sh/uv/). From this directory:

```bash
uv sync
# optional: source .venv/bin/activate
# or prefix every command with `uv run`
```

Dependencies are in `pyproject.toml`; `uv.lock` pins versions.

## Data required

Point `--data-dir` at the **06_v5_label_aln** analysis directory, or set
`ANDV_DATA_DIR`. Default if unset:

`../hantavirus/Alingments Piet/analyses/06_v5_label_aln`

(i.e. `/mnt/storage/qc2358/hantavirus/Alingments Piet/analyses/06_v5_label_aln`).

Required for `--scope full`:

- `metadata/andv_attributes_v5.csv`
- `alignments/S_v5_label_aln.fasta`
- `alignments/M_v5_label_aln.fasta`
- `alignments/L_v5_label_aln.fasta`

Required for `--scope cladeIII` (publication examples):

- `clade_III/metadata/clade_III_attributes_v5.csv`
- `clade_III/alignments/S_cladeIII_v5.fasta`
- `clade_III/alignments/M_cladeIII_v5.fasta`
- `clade_III/alignments/L_cladeIII_v5.fasta`

Do not commit large FASTA/CSV datasets into this repo. Outputs go under
`--out-dir` (defaults under `results/`). **Tracked publication figures and
tables are already in `results/`.**

## Methods (brief)

### Pairwise segment distance scatter
For every tip pair, compute observed p-distance on S, M, and L (ignoring gap/N
sites). Plot S vs M, M vs L, and S vs L. Points off the y=x diagonal have
segment-incongruent distances (possible reassortment signal).

### Sliding-window p-distance
Along comparable (non-gap/N) sites, slide a window (default 100 nt, step 10)
and plot differences / window size for S/M/L for one isolate pair.

### Mutation-position KS test
For all tip pairs, collect substitution coordinates among comparable sites and
test against **Uniform(0, n_comparable)**. Larger *D* / smaller *p* ⇒ clustered
mutations. Ranked by L-segment KS *D* (requires ≥ `--min-mut` L mutations).

## Usage

Publication examples use the **06_v5_label_aln** Clade III data. Publication
outputs (binomial scatters, sliding-window panels, KS tables)
are git-tracked under `results/`.

```bash
export ANDV_DATA_DIR="/mnt/storage/qc2358/hantavirus/Alingments Piet/analyses/06_v5_label_aln"

# Pairwise S/M/L scatters (Clade III) — underlies binomial publication panels
# Tracked: results/segment_distance_scatter/CladeIII_pairwise_{S_vs_M,S_vs_L,M_vs_L}_binomial_publication.png
uv run python scripts/plot_pairwise_segment_distance_scatter.py \
  --scope cladeIII --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/segment_distance_scatter

# Sliding window: p1236 vs p1239 (window 100, step 10)
# Tracked: results/sliding_window/p1236_vs_p1239_sliding_window_pdist_publication.*
uv run python scripts/plot_sliding_window_pdist.py \
  --seq-a p1236 --seq-b p1239 \
  --window 100 --step 10 --scope cladeIII \
  --data-dir "$ANDV_DATA_DIR" --out-dir results/sliding_window

# Sliding window: p1059 vs NRC-4/18 (window 100, step 10)
# Tracked: results/sliding_window/p1059_vs_NRC-4-18_sliding_window_pdist_publication.*
uv run python scripts/plot_sliding_window_pdist.py \
  --seq-a p1059 --seq-b "NRC-4/18" \
  --window 100 --step 10 --scope cladeIII \
  --data-dir "$ANDV_DATA_DIR" --out-dir results/sliding_window

# Mutation-position KS (all pairs + rank; publication pair in results/KS/)
# Tracked: results/KS/p1236_vs_p1239_L_mutation_positions_KS.csv
uv run python scripts/ks_mutation_positions.py \
  --scope cladeIII --min-mut 20 \
  --data-dir "$ANDV_DATA_DIR" --out-dir results/KS
```

## Parameters

### Shared

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--data-dir` | `ANDV_DATA_DIR` or `../hantavirus/Alingments Piet/analyses/06_v5_label_aln` | Root with v5 alignments and metadata |
| `--out-dir` | script-specific under `results/` | Where PNG/CSV outputs are written |
| `--scope` | see below | `full` = all v5 tips; `cladeIII` = Clade III subset |

### `plot_pairwise_segment_distance_scatter.py`

All tip pairs → S/M/L p-distances → scatter panels (no reference isolate).

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--scope` | `cladeIII` | Tip set for pairwise comparisons |
| `--highlight` | `p1236` | Substrings matching either tip of a pair → blue highlight (can pass several) |
| `--exclude` | _(none)_ | Substrings of tips to drop before computing pairs |
| `--combined-only` | off | Write only the 3-panel figure (skip per-panel S–M / M–L / S–L PNGs) |

Cruise-ship pairs are always colored red; HHPC tips use diamond markers.
Outputs: `{scope}_pairwise_segment_pdist.csv` and scatter PNG.
Tracked binomial publication panels + lambda-test CSVs live under
`results/segment_distance_scatter/`.

### `plot_sliding_window_pdist.py`

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--seq-a` | _(required)_ | First isolate (`sample_descriptor`, accession, or alias e.g. `cruise`) |
| `--seq-b` | _(required)_ | Second isolate (same resolution rules) |
| `--window` | `100` | Window length in comparable (non-gap/N) sites |
| `--step` | `10` | Step between window starts (comparable-site units) |
| `--scope` | `full` | Which alignment/metadata set to load |
| `--per-segment` | off | Also write separate S, M, L single-panel figures |

Outputs: `{A}_vs_{B}_sliding_window_pdist.{png,csv}` (+ optional `_*_{S,M,L}.*`).

### `ks_mutation_positions.py`

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--scope` | `full` | Tip set for all-pairs mutation KS |
| `--min-mut` | `20` | Minimum L mutations required for a non-blank `rank_L` |

Outputs: `mutation_positions_KS_vs_uniform_all_pairs.csv` and a short README txt.

## Repository contents

| Path | Role |
|------|------|
| `pyproject.toml` / `uv.lock` | Dependencies (uv) |
| `scripts/plot_pairwise_segment_distance_scatter.py` | All-pairs S/M/L p-distance scatters |
| `scripts/plot_sliding_window_pdist.py` | Pairwise sliding-window p-distance plots |
| `scripts/ks_mutation_positions.py` | KS of mutation positions vs Uniform |
| `results/` | Publication figures/tables (**git-tracked**) |
| `.venv/` | Local virtualenv (gitignored; `uv sync`) |
