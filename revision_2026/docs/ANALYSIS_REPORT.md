# Additional analyses for Epigenomics revision

## Dataset and scope

Existing GSE70090 annotated oxBS sample: 250,000 independently sampled CpGs from each of four patients × normal/tumor (2,000,000 observations), hg19; coverage ≥10. This is not a new experiment, independent cohort, or matched-locus normal–tumor sample. Site-level analysis is based on sampled observations.

## Leave-one-patient-out validation

The baseline predicts methylation from status × genomic compartment × coverage stratum. The augmented model adds an additive, training-only 256-NNCGNN residual with a fixed 1,000-site shrinkage. No held-out patient is used in training or tuning.

| Patient held out | Baseline RMSE (pp) | + NNCGNN RMSE (pp) | Relative reduction |
|---|---:|---:|---:|
| 1 | 30.06 | 29.52 | 1.81% |
| 2 | 33.95 | 33.60 | 1.05% |
| 3 | 30.71 | 30.15 | 1.83% |
| 4 | 32.60 | 32.20 | 1.23% |

Prediction of absolute site methylation improved modestly in all four patients. R² improved from 0.08–0.11 to 0.10–0.13 (see metrics CSV). These are patient-held-out estimates within a single four-patient cohort; no external cohort was tested.

Prediction of motif-specific tumor–normal methylation differences did **not** improve: motif-delta RMSE was 13.7109 pp with baseline and 13.7126 pp with additive context. Adding a status-specific motif residual yielded 13.7371 pp. Do not describe this as validated prediction of tumor remodeling.

## Regulatory / functional contextualization

Six mutually exclusive hg19 annotation categories are available: CpG island, non-CGI promoter, CDS, non-CDS exon, other genic/intronic, intergenic non-CGI. The existing sample contains no gene IDs, expression or chromatin-state data.

ATCGAT in non-CGI promoters: mean four-patient compartment-adjusted tumor-minus-normal residual −13.99 pp (negative in all four; approximately 30–48 tumor sites and 34–43 normal sites per patient). CGCGCG in intergenic non-CGI regions: +15.30 pp (positive in all four; 55–70 tumor sites and 52–76 normal sites per patient). Residuals represent relative retention/depletion against compartment mean, not absolute hyper- or hypomethylation. These are descriptive, selected motifs, not discovery-adjusted inferential findings.

For candidate regional localization, 1-Mb blocks were summarized separately for promoter and CpG-island strata, requiring at least 10 CpGs in each condition. 4,214 patient×block×compartment comparisons across 1,281 unique blocks passed this filter. These are coarse regions, not assigned genes or regulatory targets.

## Remaining limitations

No gene-resolved hg19 annotation table is included in the uploaded site-level sample, and no matched RNA-seq/ATAC-seq/ChIP-seq data are available here. Gene expression relationships and chromatin effects cannot be computed without obtaining those datasets. Sample size is four biological pairs; p-values based on millions of sites would overstate biological replication. Motif analyses use independently sampled normal and tumor sites and cannot be treated as matched CpG pairs. Clinical/variant-effect prediction has not been demonstrated.

## Files

- `heldout_prediction_metrics.csv`, `heldout_prediction_comparison.csv`, `heldout_compartment_metrics.csv`: full prediction benchmarking.
- `heldout_motif_tumor_delta.csv`, `heldout_status_interaction_metrics.csv`: change-prediction sensitivity.
- `functional_compartment_selected_motifs_adjusted.csv`: motif-by-patient, compartment-centered residuals.
- `regulatory_region_1Mb_patient_changes.csv`: regional exploratory results.
- `site_summary_sufficient_statistics.csv`: model input sufficient statistics.
- `run_validation.py`, `interaction_validation.py`: reproducible analysis.
- `Epigenomics_validation_exploratory_editable.pptx`: draft editable Supplementary Figure S4 panels.

To reproduce, place `GSE70090_oxBS_site_level_annotated.csv.gz` at the local source directory and run the commands in `REPRODUCE.md`. The scripts use fixed shrinkage values, no randomized model tuning, and training/test split by patient.
