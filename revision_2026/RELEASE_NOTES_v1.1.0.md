# Release notes — v1.1.0 (to publish after file upload)

October 2026 reviewer-response update for “Conserved local sequence context predicts CpG and non-CpG methylation across human methylomes.”

Adds exploratory 256-motif and genomic-compartment summaries, patient-stratified normal–tumor comparisons, leave-one-patient-out absolute methylation prediction, tumor-change prediction sensitivity, 1-Mb regional summaries, documentation and reproducible Python scripts.

The annotated input comprises 2,000,000 independently sampled CpG observations from four matched patient pairs (not matched CpG loci). Sequence context modestly improves held-out absolute methylation prediction, but not prediction of tumor–normal motif changes. This is not an independent external cohort or gene-expression validation. Large primary WGBS files are not redistributed.

**Before release:** Upload all data, scripts, and documentation from the prepared ZIP, verify their presence, merge the branch to main, update VERSION and CITATION.cff to 1.1.0, then create the v1.1.0 release. Preserve v1.0.0.
