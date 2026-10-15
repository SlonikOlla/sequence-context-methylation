# Reproduce the October 2026 revision analyses

## Inputs

Obtain GSE70090 processed oxBS data and the original hg19 annotation inputs as documented in `docs/DATASETS.md`. The scripts operate on an **existing derived input** named `GSE70090_oxBS_site_level_annotated.csv.gz`, containing 2,000,000 independently sampled site records (250,000 per patient × status). **The annotated site-level input is not included in this release**; the scripts alone do not regenerate it from GEO raw data. Consult the original repository workflows and ensure that this intermediate has been generated and verified before running these scripts.

For the flanking-position and full-motif summaries, also provide the relevant WGBS summary ZIPs (`WGBS_context_results.zip`, `WGBS_compartment_results_Liver53F.zip`, `WGBS_compartment_results_Pancreas34M.zip`, `WGBS_compartment_results_SigmoidColon.zip`, `WGBS_compartment_results_HepG2.zip`) from the earlier analysis. These ZIPs are not redistributed here.

## Python environment

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install 'pandas>=2.2,<3' 'numpy>=1.26,<3'
```

## Commands

From the repository root:

```bash
python revision_2026/scripts/run_validation.py --input /path/to/GSE70090_oxBS_site_level_annotated.csv.gz --output-dir revision_2026/data_recomputed
python revision_2026/scripts/interaction_validation.py --output-dir revision_2026/data_recomputed
python revision_2026/scripts/compute_revision_data.py --source-dir /path/to/source_archives --output-dir revision_2026/data_recomputed
```

All model coefficients are computed using training patients only. Holdout patient outcomes are never used to tune the fixed shrinkage of 1,000 sites for the pooled motif residual and 1,500 for status-specific motif sensitivity. Compare outputs against the supplied tables. Values may differ slightly by numerical precision and library versions.

## Caveats

The `site_summary_sufficient_statistics.csv` file contains the model's per-group sums and squared sums. This permits checking the held-out RMSE calculations without access to individual CpG measurements, but **not** reconstructing chromosome-level, regulatory-region or all patient-site relationships. Biological replication is four patient pairs. Do not treat the 2M site records as 2M independent patients.
