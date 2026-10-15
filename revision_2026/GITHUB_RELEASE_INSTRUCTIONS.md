# GitHub update and release checklist

Repository: https://github.com/SlonikOlla/sequence-context-methylation

1. Upload the `revision_2026/` directory to the repository root (do not upload the parent `GitHub_revision_2026/` directory).
2. Review the CSV files and script files in the new directory.
3. Update `VERSION` from `1.0.0` to `1.1.0`; update `CITATION.cff` version/date only if all authors and repository metadata are correct.
4. Commit with message `Add October 2026 reviewer-response data and validation scripts`.
5. Create a GitHub release with tag `v1.1.0` and target `main` **after** all files are present. Do not change the historical `v1.0.0` tag.

Release title: Sequence-context methylation — revised reproducibility package v1.1.0

Release notes:

October 2026 reviewer-response update for “Conserved local sequence context predicts CpG and non-CpG methylation across human methylomes.” This release adds exploratory full-motif and genomic-compartment summaries from GSE70090, patient-stratified normal–tumor comparisons, leave-one-patient-out absolute methylation prediction, a tumor-change prediction sensitivity analysis, 1-Mb regional summaries, documentation, and reproducible Python scripts. The site-level input contains 2,000,000 independently sampled CpG observations from four matched patient pairs (not matched CpG loci). Sequence context modestly improves held-out absolute methylation prediction but does not improve prediction of tumor–normal motif changes. No new independent cohort or functional gene-expression validation is claimed. Large primary WGBS files are not redistributed.
