#!/usr/bin/env python3
"""Fit adjusted grouped-binomial normal-vs-tumor sequence-context models for GSE70090.

Input: GSE70090_oxBS_site_level_annotated.csv.gz created by
annotate_GSE70090_site_samples_hg19.py.

The script aggregates sampled CpGs by patient, status, coverage band,
primary genomic compartment, and NNCGNN motif, then compares:

M0: patient + status + patient:status + coverage + compartment:status + motif
M1: M0 + status:motif

The status x motif terms test whether specific sequence contexts show relative
methylation retention or depletion in tumor after adjustment for patient,
coverage, and genomic compartment.

Requires pandas, numpy, scipy and statsmodels.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2

if len(sys.argv) < 2:
    raise SystemExit("Usage: python fit_GSE70090_cancer_context_model.py GSE70090_oxBS_site_level_annotated.csv.gz [output_dir]")

infile = Path(sys.argv[1])
outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else infile.parent / "GSE70090_cancer_model_results"
outdir.mkdir(parents=True, exist_ok=True)

usecols = ["patient","status","coverage_band","primary_compartment","sixmer","methylated_reads","coverage","block_1Mb"]
df = pd.read_csv(infile, usecols=usecols)
for c in ["patient","status","coverage_band","primary_compartment","sixmer"]:
    df[c] = df[c].astype("category")

# Aggregate read counts over identical design cells. This yields the same
# grouped-binomial likelihood while keeping the model computationally tractable.
g = (df.groupby(["patient","status","coverage_band","primary_compartment","sixmer"], observed=True)
       .agg(methylated_reads=("methylated_reads","sum"), coverage=("coverage","sum"), n_sites=("sixmer","size"))
       .reset_index())
g["meth_prop"] = g["methylated_reads"] / g["coverage"]

base = "meth_prop ~ C(patient) + C(status) + C(patient):C(status) + C(coverage_band) + C(primary_compartment)*C(status) + C(sixmer)"
full = base + " + C(status):C(sixmer)"

rows=[]
fits={}
for name,formula in [("M0_no_tumor_by_sequence",base),("M1_tumor_by_sequence",full)]:
    print("Fitting", name)
    model = smf.glm(formula=formula, data=g, family=sm.families.Binomial(), freq_weights=g["coverage"])
    fit = model.fit(maxiter=200)
    fits[name]=fit
    rows.append([name, fit.deviance, fit.aic, len(fit.params), fit.nobs])
    ci=fit.conf_int()
    tab=pd.DataFrame({"term":fit.params.index,"coef_logit":fit.params.values,"SE":fit.bse.values,"z":fit.tvalues.values,"p":fit.pvalues.values,"CI_low":ci[0].values,"CI_high":ci[1].values})
    tab.to_csv(outdir/f"{name}_coefficients.csv",index=False)

pd.DataFrame(rows,columns=["model","deviance","AIC","n_parameters","n_design_cells"]).to_csv(outdir/"model_comparison.csv",index=False)
m0,m1=fits["M0_no_tumor_by_sequence"],fits["M1_tumor_by_sequence"]
dev=m0.deviance-m1.deviance
df_diff=len(m1.params)-len(m0.params)
p=chi2.sf(max(dev,0),df_diff)
pd.DataFrame([[dev,df_diff,p]],columns=["deviance_reduction","df","LR_p"]).to_csv(outdir/"tumor_by_sequence_deviance_test.csv",index=False)

# Export just tumor-by-sixmer coefficients, sorted by effect.
coef=pd.read_csv(outdir/"M1_tumor_by_sequence_coefficients.csv")
mask=coef.term.str.contains("C\(status\).*:C\(sixmer\)|C\(sixmer\).*:C\(status\)", regex=True)
seq=coef[mask].sort_values("coef_logit")
seq.to_csv(outdir/"tumor_by_sixmer_coefficients_sorted.csv",index=False)

print("Done:", outdir)
