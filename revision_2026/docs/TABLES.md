# Supplementary Data — table inventory

Manuscript: Conserved local sequence context predicts CpG and non-CpG methylation across human methylomes.
Revision: October 2026; derived data package 1.1.0.

| File | Data rows | Key columns |
|---|---:|---|
| `cancer_motif_compartment_four_patient_summary.csv` | 1,019 | primary_compartment, sixmer, n_patients, mean_residual, sd_residual |
| `cancer_motif_compartment_patient_effects.csv` | 6,144 | patient, primary_compartment, sixmer, n_Normal, n_Tumor |
| `cancer_motif_four_patient_global_adjusted.csv` | 1,024 | patient, sixmer, n_Normal, n_Tumor, mean_Normal |
| `full_flanking_position_effects.csv` | 120 | dataset, position, base, mean_methylation_percent, difference_from_A_pp |
| `functional_compartment_selected_motifs.csv` | 144 | patient, primary_compartment, sixmer, n_Normal, n_Tumor |
| `functional_compartment_selected_motifs_adjusted.csv` | 144 | patient, primary_compartment, sixmer, n_Normal, n_Tumor |
| `heldout_compartment_metrics.csv` | 48 | patient, compartment, model, n, RMSE_pp |
| `heldout_motif_tumor_delta.csv` | 1,024 | patient, sixmer, n_Normal, n_Tumor, observed_Normal |
| `heldout_prediction_comparison.csv` | 12 | held_out_patient, subset, baseline, plus_6mer, RMSE_improvement_pp |
| `heldout_prediction_metrics.csv` | 24 | held_out_patient, subset, model, n_sites, RMSE_pp |
| `heldout_status_interaction_metrics.csv` | 12 | patient, subset, RMSE_pp |
| `regulatory_region_1Mb_patient_changes.csv` | 4,214 | patient, block_1Mb, primary_compartment, n_Normal, n_Tumor |
| `site_summary_sufficient_statistics.csv` | 48,317 | patient, status, primary_compartment, coverage_band, group |
| `train_only_motif_effects.csv` | 1,024 | held_out_patient, sixmer, n, effect |

## Definitions and interpretation

- `sixmer`: 6-bp NNCGNN sequence context; 256 possible combinations of four flanking positions.
- `patient`: one of four biological normal–tumor liver pairs; a patient, not a CpG, is the unit of biological replication.
- `status`: Normal or Tumor.
- `primary_compartment`: mutually exclusive hg19 genomic compartment defined in original annotation.
- `coverage_band`: sequencing-depth category used in the fitted baseline.
- `n`, `n_Normal`, `n_Tumor`: sampled site counts, not numbers of patients.
- `mean`, `observed`, `sy`: methylation fraction on 0–1 scale unless column explicitly ends in `_pp` or `_percent`.
- `sy2`: sum of squared methylation fractions for RMSE computation.
- `delta_pp`: tumor minus normal mean methylation in percentage points.
- `residual_delta_pp`: motif tumor–normal change minus the same patient and compartment mean change, in percentage points.
- `RMSE_pp`: root mean squared prediction error in methylation percentage points; lower is better.
- `relative_RMSE_reduction_pct`: percentage decrease in RMSE vs baseline; this is not percentage-point change.
- `block_1Mb`: chromosome-specific 1-Mb block label; not a gene assignment.
- `difference_from_A_pp`: descriptive mean methylation difference from base A at the same flank position; not a causal allele substitution estimate.

## Scope and limitations

- Source: GSE70090, hg19/GRCh37, four matched normal–tumor liver patients; the annotated reservoir contains 250,000 oxBS site observations per patient × condition (2,000,000 total), filtered at coverage ≥10.
- Normal and tumor observations are **not matched at identical CpG coordinates**; the cohort is paired by patient only.
- The validation is leave-one-patient-out **within the same cohort**; no external patient cohort is claimed.
- The additive 256-motif model modestly improves prediction of absolute site methylation in all four held-out patients, but does **not** improve prediction of motif-specific tumor–normal changes.
- Selected motifs and 1-Mb blocks are exploratory, not inferentially adjusted discoveries or demonstrated gene regulatory mechanisms.
- No matched RNA expression, chromatin or gene ID assignments are provided by these tables.
- Primary WGBS/oxBS data and genome reference files are excluded. Consult GEO and the repository dataset manifest.

## Reproducibility

Analysis scripts and detailed methods are in `revision_2026/scripts/` in the companion GitHub package. Scripts require Python 3, pandas and numpy; see `revision_2026/docs/REPRODUCE.md`.

## Citation

Cite the manuscript and the versioned GitHub repository release when available. Do not cite this archive as a new independent cohort.
