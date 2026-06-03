# HAWK-EYE exon annotation pipeline

This repository contains the scripts and outputs used in the HAWK-EYE exon-level analysis. The included Snakemake workflow (`Snakefile`) provides a reproducible run order that mirrors the manual steps used for the published analysis.

Top-level steps (as implemented in the `Snakefile`):

- `1-Parse ExonCalculator`: parse ExonCalculator outputs with `parse_exoncalculator.R`.
- `2-ClinVar`: annotate ClinVar LPP pathogenic variants with `lpp_clinvar_annotation.R` (requires ClinVar files in `2-ClinVar/clinvar_db`).
- `3-Premature Stop`: compute exon sequence features using CCDS with `exon_features.R`.
- `4-PhyloP`: compute exon-level phyloP scores with `phylop_mean_conservation_score.py` (requires a phyloP bigWig in `4-PhyloP/phylop/`).
- `5-UniProt` → `6-InterPro` → `7-GTEx` → `8-PSI` → `9-Domain Counting`: downstream annotation steps implemented as R scripts.

Prerequisites
- Install Snakemake (recommended via conda or mamba):

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba -y
conda activate snakemake
```

- Ensure `Rscript` and `python3` are available. To use per-rule conda environments, Snakemake will create the envs from the `envs/` YAML files.

Running the full workflow

From the repository root run:

```bash
snakemake --use-conda -j 4
```

This will create per-rule conda environments under `envs/` as specified, and run the pipeline in the correct order.

Running a single rule (example)

```bash
snakemake '4-PhyloP/HAWKEYE_Database_PhyloP_Exon_Conservation.csv' --use-conda -j 2
```

Data availability notes
- Large or restricted input files (ClinVar dumps, CCDS FASTA files, and the phyloP bigWig) are not included in this repository. Place them in the paths expected by the scripts: see `config.yaml` and the per-folder `clinvar_db/` and `phylop/` subdirectories.

Next steps / recommendations
- Validate the `Snakefile` on a small example dataset (or a subset of your data) using `-j 1`.
- Consider adding a small example input set under `tests/data/` for CI.
- Optionally add a GitHub Actions workflow that runs a tiny end-to-end example.

If you want, I can:

- refine the `Snakefile` to capture per-file parameterization and add per-rule logs,
- add a small `tests/` example and a GitHub Actions CI workflow,
- or create a `Dockerfile` for runtime encapsulation.
