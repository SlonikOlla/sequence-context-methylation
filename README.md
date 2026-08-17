# Sequence-context analysis of human DNA methylation

Reproducibility package for the manuscript **“Conserved local sequence context predicts CpG and non-CpG methylation across human methylomes.”**

This repository contains the analysis code, workflow documentation, public-dataset manifest, reduced summary tables, and editable draft figures needed to reproduce the analyses without redistributing the very large primary WGBS files.

## Contents

- `scripts/` – Python scripts for CpG/CHG sequence-context reduction, coverage stratification, genomic annotation, motor-neuron CpH analysis, matched normal–tumor analysis, and multivariable modeling.
- `docs/DATASETS.md` – public accession manifest and required reference annotations.
- `docs/WORKFLOW.md` – command-line workflow and analysis order.
- `docs/GITHUB_UPLOAD_INSTRUCTIONS.md` – exact GitHub repository text and upload instructions.
- `reduced_results/` – compact result tables used in the manuscript analyses. Large primary WGBS files are excluded.
- `figures/` – editable PowerPoint versions of the current figures.
- `requirements.txt` – Python dependencies used for statistical modeling and figure-generation utilities.
- `CITATION.cff` – citation metadata recognized by GitHub.
- `LICENSE` – MIT license for the code in this repository.
- `SHA256SUMS.txt` – checksums for repository files in this release package.

## Reproducibility model

Primary sequencing/methylation files are downloaded from their public repositories using the accessions listed in `docs/DATASETS.md`. The reduction scripts convert those large files into compact sequence-context tables. Subsequent scripts perform coverage-stratified analyses, genomic-compartment annotation, neuronal CpH validation, normal–tumor comparisons, and adjusted modeling.

The repository intentionally does **not** redistribute large ENCODE/GEO WGBS files, reference genomes, or annotation files.

## Python

The core reduction scripts use the Python standard library. Statistical modeling uses pandas, NumPy, SciPy, statsmodels and patsy. See `requirements.txt`.

## Before making the GitHub repository public

1. Replace `REPLACE_WITH_GITHUB_USERNAME` in `CITATION.cff` with the GitHub account or organization name that will host the repository.
2. Add all manuscript authors, in final author order, to `CITATION.cff`; add ORCID identifiers where available.
3. Confirm the manuscript title in `README.md` and `CITATION.cff` after the final proofread.
4. If desired, add local acquisition dates and SHA256 checksums for the large downloaded source files to `docs/DATASETS.md`. The public accessions themselves are already recorded.
5. Upload the **contents of this folder** to the root of the GitHub repository.

## Citation

Please cite the associated manuscript. GitHub will also display citation information from `CITATION.cff` once the repository is public.
