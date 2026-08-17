#!/usr/bin/env python3
"""Fit site-level binomial GLMs to output from prepare_site_level_cpg.py.

Requires: pandas, numpy, scipy, statsmodels, patsy.
Model 0: sample + coverage + genomic compartment
Model 1: Model 0 + six individual flanking bases
Model 2: Model 0 + complete NNCGNN sixmer
Also fits a focused interaction model testing sample-specific -1/+1/+3 effects.
Cluster-robust standard errors use 1-Mb genomic blocks.
"""
import sys, pandas as pd, numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2
from pathlib import Path

if len(sys.argv)<2:
    raise SystemExit("Usage: python fit_site_level_cpg_model.py cpg_site_level_sample.tsv.gz [output_dir]")
infile=Path(sys.argv[1]);out=Path(sys.argv[2]) if len(sys.argv)>2 else infile.parent/"model_results";out.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(infile,sep="\t")
df["meth_prop"]=df.methylated_reads/(df.methylated_reads+df.unmethylated_reads)
df["log_cov"]=np.log1p(df.coverage)
# Sampling weights restore stratum representation; coverage is separately modeled, not used as an extra frequency weight.
w=df.sampling_weight.astype(float)
formulas={
"M0":"meth_prop ~ C(sample) + log_cov + C(compartment)",
"M1_positions":"meth_prop ~ C(sample) + log_cov + C(compartment) + C(m3)+C(m2)+C(m1)+C(p1)+C(p2)+C(p3)",
"M2_sixmer":"meth_prop ~ C(sample) + log_cov + C(compartment) + C(sixmer)",
"M3_focused_interactions":"meth_prop ~ C(sample) + log_cov + C(compartment) + C(m3)+C(m2)+C(m1)+C(p1)+C(p2)+C(p3) + C(sample):C(m1) + C(sample):C(p1) + C(sample):C(p3)"
}
fits={};rows=[]
for name,formula in formulas.items():
    print("Fitting",name)
    model=smf.glm(formula=formula,data=df,family=sm.families.Binomial(),freq_weights=w)
    fit=model.fit(cov_type="cluster",cov_kwds={"groups":df["block_1Mb"]},maxiter=100)
    fits[name]=fit
    rows.append([name,fit.nobs,fit.llf,fit.deviance,fit.aic,len(fit.params)])
    ci=fit.conf_int();tab=pd.DataFrame({"term":fit.params.index,"coef_logit":fit.params.values,"SE_cluster":fit.bse.values,"z":fit.tvalues.values,"p":fit.pvalues.values,"CI_low":ci[0].values,"CI_high":ci[1].values})
    tab.to_csv(out/f"{name}_coefficients.csv",index=False)
pd.DataFrame(rows,columns=["model","n_rows","logLik","deviance","AIC","n_parameters"]).to_csv(out/"model_comparison.csv",index=False)
# Deviance improvements relative to M0 (descriptive LRT; robust SE are in coefficient tables).
base=fits["M0"]
lr=[]
for name in ["M1_positions","M2_sixmer"]:
    f=fits[name];stat=base.deviance-f.deviance;dfdiff=len(f.params)-len(base.params);p=chi2.sf(max(stat,0),dfdiff)
    lr.append([name,stat,dfdiff,p])
pd.DataFrame(lr,columns=["model_vs_M0","deviance_reduction","df","LR_p"]).to_csv(out/"deviance_tests.csv",index=False)
print("Done:",out)
