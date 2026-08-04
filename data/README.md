# Input data

This directory is a placeholder. Scripts expect `--data-dir` to point at the
parent **Alingments Piet** directory (typo retained from the source project),
or set `ANDV_DATA_DIR`.

Default if unset: `../Alingments Piet` relative to this package root.

## Expected input files

Required:

- `andv_attributes_joined_metadata.csv`
- `S_complete_aln_CURATED.fasta`
- `M_complete_aln_CURATED.fasta`
- `L_complete_aln_CURATED.fasta`

Optional (ML-distance scatters / Clade III scope):

- `phylogeny/iqtree/{S,M,L}/{S,M,L}_ML.mldist`
- `clade_III/metadata/clade_III_attributes_joined.csv`
- `clade_III/alignments/{S,M,L}_cladeIII_subset_from_CURATED.fasta`
- `clade_III/iqtree/{S,M,L}/{S,M,L}_ML.mldist`

Do not commit large FASTA/CSV datasets into this repo (see root `.gitignore`).
