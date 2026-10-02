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

Point `--data-dir` (or `ANDV_DATA_DIR`) at a directory with this layout. Filenames are fixed.

`--scope full`:

- `metadata/andv_attributes_v5.csv`
- `alignments/S_v5_label_aln.fasta`
- `alignments/M_v5_label_aln.fasta`
- `alignments/L_v5_label_aln.fasta`

`--scope cladeIII`:

- `clade_III/metadata/clade_III_attributes_v5.csv`
- `clade_III/alignments/S_cladeIII_v5.fasta`
- `clade_III/alignments/M_cladeIII_v5.fasta`
- `clade_III/alignments/L_cladeIII_v5.fasta`

### Metadata CSV

Comma-separated, one row per tip. Both scopes use the same columns. Extra columns are ignored.

| Column | Required | Role |
|--------|----------|------|
| `sample_descriptor` | yes | Tip name used in outputs and in `--seq-a`, `--seq-b`, `--exclude`, `--fit-exclude`, and `--plot-exclude` |
| `full_label` | yes | Alternate name for sample lookup and short labels (for example Chile-9717869) |
| `meta_S`, `meta_M`, `meta_L` | yes | FASTA record id for that tip on each segment. Every value must exist in the matching alignment |
| `sequence_group` | no | A value starting with `Cruise` is the `cruise` alias |
| `andv_clade` | no | Clade label (for example `III`). Used to mark Clade III pairs on full-scope scatters and in KS output |
| `highlighted_HHPC_cluster` | no | Copied to the pairwise table as `hhpc_a` / `hhpc_b` |

### Alignments

Three FASTA files, one per segment. Within a file every sequence has the same length (a multiple alignment). S, M, and L may have different lengths.

The record id is the first token of the header, up to the first whitespace. That id must equal `meta_S`, `meta_M`, or `meta_L` for the corresponding segment. A header `>PZ485159_127D` has id `PZ485159_127D`.

Comparable sites are columns where both sequences are `A`, `C`, `G`, or `T` (case is ignored). `-` and `N` are skipped for p-distance, sliding windows, and the KS test. Any other character is treated as a called base.

## Methods

### Pairwise segment distance scatter (binomial λ test)
For every tip pair, compute mutation counts and comparable-site opportunities
on S, M, and L. Fit a conditional binomial rate-ratio λ
(*Y<sub>Y</sub>* | *T* ~ Binomial(*T*, π), π = λ *E<sub>Y</sub>* / (*E<sub>X</sub>* + λ *E<sub>Y</sub>*))
on pairs with expected count *T*π > `--mu-min`. Test each pair with a normal
approximation Z-test; apply Benjamini–Hochberg FDR at `--fdr-alpha`.
Publication panels show *y* = λ*x*. On the full-tip panels, blue marks a pair in which both tips are Clade III, and a black edge marks BH *q* < 0.05.

### Sliding-window p-distance
Along comparable (non-gap/N) sites, slide a window (default 100 nt, step 10)
and plot differences / window size for S/M/L for one isolate pair.

### Mutation-position KS test
For all tip pairs, collect substitution coordinates among comparable sites and
test against **Uniform(0, n_comparable)**. Larger *D* / smaller *p* ⇒ clustered
mutations. Ranked by L-segment KS *D* (requires ≥ `--min-mut` L mutations).

## Usage

Retained results use `07_L_v5_label_aln_3_no_R600` (new L alignment, R600 excluded).

```bash
export ANDV_DATA_DIR="/mnt/storage/qc2358/hantavirus/Alignment/analyses/07_L_v5_label_aln_3_no_R600"

# Clade III pairwise scatters
uv run python scripts/plot_pairwise_segment_distance_scatter.py \
  --scope cladeIII --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/segment_distance_scatter_07_no_R600 \
  --mu-min 20 --fdr-alpha 0.05

# All tips, R600 excluded
uv run python scripts/plot_pairwise_segment_distance_scatter.py \
  --scope full --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/segment_distance_scatter_full_Laln3_no_R600 \
  --mu-min 20 --fdr-alpha 0.05

# Sliding windows
uv run python scripts/plot_sliding_window_pdist.py \
  --seq-a p1236 --seq-b p1239 \
  --window 100 --step 10 --scope cladeIII --publication \
  --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/sliding_window_L_v5_label_aln_3

uv run python scripts/plot_sliding_window_pdist.py \
  --seq-a p1059 --seq-b "NRC-4/18" \
  --window 100 --step 10 --scope cladeIII --publication \
  --data-dir "$ANDV_DATA_DIR" \
  --out-dir results/sliding_window_L_v5_label_aln_3
```

The one-pair L KS table for p1236 vs p1239 is `results/KS_L_v5_label_aln_3/p1236_vs_p1239_L_mutation_positions_KS.csv`.


## Parameters

### Shared

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--data-dir` | `ANDV_DATA_DIR` or `../hantavirus/Alignment/analyses/06_v5_label_aln` | Root with v5 alignments and metadata |
| `--out-dir` | script-specific under `results/` | Where PNG/CSV outputs are written |
| `--scope` | see below | `full` = all v5 tips; `cladeIII` = Clade III subset |

### `plot_pairwise_segment_distance_scatter.py`

All tip pairs → mutation counts / opportunities → conditional binomial λ test →
publication scatters.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--scope` | `cladeIII` | Tip set for pairwise comparisons |
| `--exclude` | _(none)_ | Tips to drop before any pairwise comparisons |
| `--fit-exclude` | _(none)_ | Tips omitted from the λ fit but kept in the BH test and, unless also plot-excluded, in the scatter |
| `--plot-exclude` | _(none)_ | Tips omitted from the scatter only. They remain in the BH test |
| `--mu-min` | `20` | Minimum expected Y count *T*π for λ fit and BH testing |
| `--fdr-alpha` | `0.05` | Benjamini–Hochberg FDR threshold (significant if *q* < alpha) |

Outputs: `{prefix}_pairwise_segment_pdist.csv`,
`{prefix}_{X}_vs_{Y}_lambda_binomial_test.csv`, and
`{prefix}_pairwise_{X}_vs_{Y}_binomial_publication.png` and `.svg`.

### `plot_sliding_window_pdist.py`

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--seq-a` | _(required)_ | First isolate (`sample_descriptor`, accession, or alias e.g. `cruise`) |
| `--seq-b` | _(required)_ | Second isolate (same resolution rules) |
| `--window` | `100` | Window length in comparable (non-gap/N) sites |
| `--step` | `10` | Step between window starts (comparable-site units) |
| `--scope` | `full` | Which alignment/metadata set to load |
| `--per-segment` | off | Also write separate S, M, L single-panel figures |
| `--publication` | off | Stacked publication panel (`*_publication.png` and `.svg`) |

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
| `results/segment_distance_scatter_07_no_R600/` | Clade III binomial scatters, new L, no R600 |
| `results/segment_distance_scatter_full_Laln3_no_R600/` | All-tip binomial scatters, new L, no R600 |
| `results/sliding_window_L_v5_label_aln_3/` | Publication sliding windows for p1236 vs p1239 and p1059 vs NRC-4/18 |
| `results/KS_L_v5_label_aln_3/` | p1236 vs p1239 L-segment mutation-position KS |
