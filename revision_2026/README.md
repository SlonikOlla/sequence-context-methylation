# October 2026 revision — additional analyses (v1.1.0 candidate)

This directory contains derived results and scripts for the revised manuscript and the responses to reviewers. It supplements, rather than replaces, the original v1.0.0 analyses. The original GitHub release v1.0.0 must remain unchanged.

- `data/`: 14 derived CSV tables, including all-motif patient/compartment summaries and held-out-patient validation.
- `scripts/`: portable command-line wrappers for the three analysis scripts.
- `docs/TABLES.md`: per-table inventory, column definitions, and scientific caveats.
- `docs/ANALYSIS_REPORT.md`: full validation findings and limitations.
- `docs/REPRODUCE.md`: required inputs and exact commands.

**Important:** no independent external cohort, matched CpG-locus comparison, gene-expression integration or validated tumor-remodeling prediction is claimed. The raw annotated 2M-site input is not redistributed and must be regenerated or obtained separately to rerun the scripts.

To make the revised submission citable, upload this entire `revision_2026/` directory, verify the tables, update the repository version to 1.1.0 and publish a GitHub `v1.1.0` release. Do not overwrite `v1.0.0`.
