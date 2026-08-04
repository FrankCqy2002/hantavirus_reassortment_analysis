# Input data

This directory is a placeholder. Scripts expect `--data-dir` to point at the
**06_v5_label_aln** analysis directory, or set `ANDV_DATA_DIR`.

Default if unset:

`../hantavirus/Alingments Piet/analyses/06_v5_label_aln`

relative to this package root (i.e.
`/mnt/storage/qc2358/hantavirus/Alingments Piet/analyses/06_v5_label_aln`).

## Expected input files

Required (`--scope full`):

- `metadata/andv_attributes_v5.csv`
- `alignments/S_v5_label_aln.fasta`
- `alignments/M_v5_label_aln.fasta`
- `alignments/L_v5_label_aln.fasta`

Required (`--scope cladeIII`):

- `clade_III/metadata/clade_III_attributes_v5.csv`
- `clade_III/alignments/{S,M,L}_cladeIII_v5.fasta`

Do not commit large FASTA/CSV datasets into this repo (see root `.gitignore`).
Publication outputs are tracked under `../results/`.
