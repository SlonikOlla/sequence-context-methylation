# October 2026 reviewer-response analyses — v1.1.0 candidate

This update supplements the original v1.0.0 release associated with the first manuscript submission. **Do not overwrite or retag v1.0.0.**

## Added analyses

- All 256 NNCGNN contexts and cytosine-flanking position summaries.
- Patient-stratified, compartment-adjusted normal–tumor liver methylation summaries.
- Leave-one-patient-out validation of absolute CpG methylation prediction, with a baseline incorporating tumor status, genomic compartment and coverage stratum.
- Sensitivity analysis of motif-specific tumor–normal change prediction.
- Exploratory 1-Mb regulatory-region summaries.

## Study limitations

The liver cohort comprises **four matched patient pairs**. The 2,000,000 annotated oxBS site observations were sampled separately from each patient × condition (250,000 per library) and **are not matched normal–tumor CpG coordinates**. Validation holds out one patient at a time within the same cohort, not an independent external cohort. Addition of sequence context modestly improves absolute methylation prediction for each held-out patient but **does not improve motif-specific tumor–normal change prediction**. Selected motifs and regions are exploratory. No gene-expression or chromatin validation is claimed.

## Files to upload before release

The companion `Epigenomics_GitHub_v1.1.0_upload.zip` contains `revision_2026/data/` (14 CSVs), `revision_2026/scripts/` (three Python scripts), and `revision_2026/docs/` (table definitions, analysis report, and reproduction instructions). Upload these files to this branch, review, and merge to main before tagging v1.1.0.

The annotated site-level input and primary WGBS/oxBS files are **not redistributed**. See `docs/DATASETS.md` for original source accessions.
