# Command-line workflow

The commands below show the intended Windows command-prompt workflow. Paths are examples and should be replaced with the local paths to the public source files.

## 1. Genome motif scan

```text
python scan_grch38_context.py
```

## 2. CpG context and genomic-compartment analyses

Run once per CpG bedMethyl file:

```text
python annotate_wgbs_compartments_v2.py
python coverage_stratify_cpg.py
```

Inputs: GRCh38 FASTA, ENCODE CpG bedMethyl, GENCODE v50 GTF, UCSC hg38 CpG-island table.

## 3. H1 CHG analyses

Run once for each H1 CHG replicate:

```text
python reduce_wgbs_chg_context.py
python coverage_stratify_chg.py
```

## 4. Motor-neuron CpH validation

Place the four ENCODE BED files in one directory and run:

```text
python reduce_motor_neuron_cph_context.py
```

Expected files:
`ENCFF093TBN.bed.gz`, `ENCFF109XWQ.bed.gz`, `ENCFF505WSM.bed.gz`, `ENCFF437XEW.bed.gz`.

## 5. Matched normal-tumor liver analysis

With `hg19.fa.gz` and `GSE70090_RAW.tar`:

```text
python reduce_GSE70090_liver_normal_tumor.py
```

Then annotate the 250,000-site reservoir samples from each library:

```text
python annotate_GSE70090_site_samples_hg19.py
```

Inputs: GENCODE v19 GTF, UCSC hg19 CpG-island table, and the `GSE70090_liver_context_results` folder.

Fit the adjusted tumor-by-sequence model:

```text
python fit_GSE70090_cancer_context_model.py GSE70090_oxBS_site_level_annotated.csv.gz
```

## 6. Cross-dataset site-level CpG model

Prepare a balanced CpG site-level sample:

```text
python prepare_site_level_cpg.py
```

Fit nested models:

```text
python fit_site_level_cpg_model.py cpg_site_level_sample.tsv.gz
```

## 7. Figures and tables

The `reduced_results/` directory contains the compact tables used for the draft figures. Editable PowerPoint versions of the current figures are in `figures/`. Before final submission, the plotting code used for the final journal-layout figures should be frozen in this repository together with the final reduced tables.
