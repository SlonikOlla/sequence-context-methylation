"""Independent patient holdout benchmarking and functional-compartment characterization.
Inputs: previously annotated reservoir samples (250k CpGs per patient x condition).
No new sequencing data or external gene-expression data.
"""
import pandas as pd, numpy as np, os
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description="Leave-one-patient-out sequence context validation")
parser.add_argument("--input", required=True, help="Annotated 2M-row oxBS site sample CSV.gz")
parser.add_argument("--output-dir", default=".", help="Directory for CSV outputs")
args=parser.parse_args()
P=Path(args.output_dir).resolve();P.mkdir(parents=True,exist_ok=True)
SRC=args.input
cols=['patient','status','chrom','position0','coverage','coverage_band','sixmer','methylation_percent','primary_compartment','CpG_island','promoter','genic','block_1Mb']
chunks=[]
for c in pd.read_csv(SRC,usecols=cols,chunksize=250000):
 c['y']=c.methylation_percent.astype('float32')/100
 c['y2']=c.y*c.y
 c['group']=c.status+'|'+c.primary_compartment+'|'+c.coverage_band
 chunks.append(c)
d=pd.concat(chunks,ignore_index=True)
print('Read',len(d),'patients',d.patient.value_counts().to_dict(), 'groups',d.group.nunique(),flush=True)
# exact sum statistics at group x motif level for reproducible out-of-sample site-level metrics
agg=d.groupby(['patient','status','primary_compartment','coverage_band','group','sixmer'],observed=True).agg(n=('y','size'),sy=('y','sum'),sy2=('y2','sum')).reset_index()
agg.to_csv(P/'site_summary_sufficient_statistics.csv',index=False)
# Bayesian-shrinkage residual motif effect learned on training patients only.
# Base: status x compartment x coverage empirical mean; augmented: base + centered 6-mer residual.
results=[]; effects=[]; comp=[]; motifdelta=[]
for hold in sorted(agg.patient.unique()):
 train=agg[agg.patient!=hold].copy(); test=agg[agg.patient==hold].copy()
 base=train.groupby('group').agg(n=('n','sum'),sy=('sy','sum')); base['base']=base.sy/base.n
 train['base']=train.group.map(base.base)
 train['residual_sum']=train.sy-train.n*train.base
 motif=train.groupby('sixmer').agg(n=('n','sum'),sr=('residual_sum','sum'))
 # shrinkage pseudocount 1000 training CpGs for each motif, not tuned on test set
 motif['effect']=motif.sr/(motif.n+1000)
 test['base']=test.group.map(base.base)
 test['effect']=test.sixmer.map(motif.effect).fillna(0)
 test['pred_base']=test.base.clip(0,1)
 test['pred_context']=(test.base+test.effect).clip(0,1)
 for k,pred in [('baseline','pred_base'),('plus_6mer','pred_context')]:
  for label,t in [('all',test),('Normal',test[test.status=='Normal']),('Tumor',test[test.status=='Tumor'])]:
   n=t.n.sum(); sq=(t.sy2-2*t[pred]*t.sy+t.n*t[pred]**2).sum();
   mae=None
   mse=sq/n
   # R2 against held-out patient global site mean
   ybar=t.sy.sum()/n
   sst=(t.sy2-2*ybar*t.sy+t.n*ybar**2).sum()
   results.append({'held_out_patient':int(hold),'subset':label,'model':k,'n_sites':int(n),'RMSE_pp':100*np.sqrt(mse),'R2':1-sq/sst,'MSE':mse})
 effects.append(motif.assign(held_out_patient=hold).reset_index()[['held_out_patient','sixmer','n','effect']])
 # independently evaluated observed motif mean per status and predicted means (aggregated over site distributions)
 for status in ['Normal','Tumor']:
  t=test[test.status==status]
  z=t.groupby('sixmer').apply(lambda q:pd.Series({'n':q.n.sum(),'observed':q.sy.sum()/q.n.sum(),'baseline':(q.n*q.pred_base).sum()/q.n.sum(),'context':(q.n*q.pred_context).sum()/q.n.sum()}),include_groups=False).reset_index()
  z['status']=status;z['patient']=hold;motifdelta.append(z)
 for compartment,t in test.groupby('primary_compartment'):
  for k,pred in [('baseline','pred_base'),('plus_6mer','pred_context')]:
   n=t.n.sum();sq=(t.sy2-2*t[pred]*t.sy+t.n*t[pred]**2).sum()
   comp.append({'patient':hold,'compartment':compartment,'model':k,'n':n,'RMSE_pp':100*np.sqrt(sq/n)})
res=pd.DataFrame(results);res.to_csv(P/'heldout_prediction_metrics.csv',index=False)
pd.concat(effects).to_csv(P/'train_only_motif_effects.csv',index=False)
pd.DataFrame(comp).to_csv(P/'heldout_compartment_metrics.csv',index=False)
m=pd.concat(motifdelta);md=m.pivot(index=['patient','sixmer'],columns='status',values=['n','observed','baseline','context'])
md.columns=['_'.join(t) for t in md.columns];md=md.reset_index()
for k in ['observed','baseline','context']:md[k+'_delta_pp']=100*(md[k+'_Tumor']-md[k+'_Normal'])
md.to_csv(P/'heldout_motif_tumor_delta.csv',index=False)
# summary of out-of-sample improvement
piv=res.pivot(index=['held_out_patient','subset'],columns='model',values='RMSE_pp').reset_index()
piv['RMSE_improvement_pp']=piv.baseline-piv.plus_6mer
piv['relative_RMSE_reduction_pct']=100*piv.RMSE_improvement_pp/piv.baseline
piv.to_csv(P/'heldout_prediction_comparison.csv',index=False)
# Regulatory / functional compartment context: patient-stratified difference for selected motifs
selected=['CGCGCG','ATCGAT','TACGTA','TCCGTA','TACGGA','AACGTT']
q=d[d.sixmer.isin(selected)].groupby(['patient','status','primary_compartment','sixmer']).agg(n=('y','size'),mean=('y','mean')).reset_index()
q=q.pivot(index=['patient','primary_compartment','sixmer'],columns='status',values=['n','mean']);q.columns=['_'.join(x) for x in q.columns];q=q.reset_index()
q['tumor_minus_normal_pp']=100*(q.mean_Tumor-q.mean_Normal)
q.to_csv(P/'functional_compartment_selected_motifs.csv',index=False)
# patient compartment baseline for residuals
b=d.groupby(['patient','status','primary_compartment']).agg(n=('y','size'),mean=('y','mean')).reset_index()
b=b.pivot(index=['patient','primary_compartment'],columns='status',values=['n','mean']);b.columns=['_'.join(x) for x in b.columns];b=b.reset_index();b['compartment_delta_pp']=100*(b.mean_Tumor-b.mean_Normal)
q=q.merge(b[['patient','primary_compartment','compartment_delta_pp']],on=['patient','primary_compartment']);q['residual_delta_pp']=q.tumor_minus_normal_pp-q.compartment_delta_pp
q.to_csv(P/'functional_compartment_selected_motifs_adjusted.csv',index=False)
# 1Mb regions: all-motif tumor-normal change for regulatory compartments with minimum sites per arm per patient
reg=d[d.primary_compartment.isin(['promoter_non_CGI','CpG_island'])]
r=reg.groupby(['patient','status','block_1Mb','primary_compartment']).agg(n=('y','size'),sy=('y','sum')).reset_index()
r=r.pivot(index=['patient','block_1Mb','primary_compartment'],columns='status',values=['n','sy']);r.columns=['_'.join(x) for x in r.columns];r=r.reset_index()
r=r.dropna();r=r[(r.n_Normal>=10)&(r.n_Tumor>=10)].copy()
r['delta_pp']=100*(r.sy_Tumor/r.n_Tumor-r.sy_Normal/r.n_Normal)
r.to_csv(P/'regulatory_region_1Mb_patient_changes.csv',index=False)
print('\nHeld-out RMSE (percentage points):\n',piv.to_string(index=False),flush=True)
print('\nRegulatory blocks with >=10 sites each condition:',len(r),'unique blocks',r.block_1Mb.nunique(),flush=True)
print('\nFunctional motifs by compartment (mean across 4 patients):')
print(q.groupby(['primary_compartment','sixmer']).residual_delta_pp.mean().unstack().round(2).to_string(),flush=True)
