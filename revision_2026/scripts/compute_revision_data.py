import pandas as pd,numpy as np,zipfile,os
import argparse
parser=argparse.ArgumentParser(description='Compute exploratory motif and compartment summaries')
parser.add_argument('--source-dir', required=True, help='Folder containing WGBS context ZIPs and annotated oxBS CSV.gz')
parser.add_argument('--output-dir', default='.', help='Directory for CSV outputs')
args=parser.parse_args()
from pathlib import Path
SOURCE_DIR=Path(args.source_dir).resolve()
O=Path(args.output_dir).resolve();O.mkdir(parents=True,exist_ok=True)
# genome-wide position contrast (base compared to A) for five datasets
mapping={'H1':'WGBS_context_results.zip','Liver':'WGBS_compartment_results_Liver53F.zip','Pancreas':'WGBS_compartment_results_Pancreas34M.zip','Sigmoid colon':'WGBS_compartment_results_SigmoidColon.zip','HepG2':'WGBS_compartment_results_HepG2.zip'}
rows=[]
for name,zipname in mapping.items():
 z=zipfile.ZipFile(SOURCE_DIR/zipname)
 target='position_summary.csv' if name=='H1' else 'compartment_position.csv'
 fn=next(n for n in z.namelist() if n.endswith(target))
 df=pd.read_csv(z.open(fn))
 print(name,df.columns.tolist(),df.head(1).to_dict('records'))
 if 'category' in df.columns:
  df=df[df.category=='ALL']
 for pos in [-3,-2,-1,1,2,3]:
  ss=df[df.position.astype(str).str.replace('+','',regex=False)==str(pos)]
  if len(ss)==0: continue
  a=ss[ss.base=='A']
  if len(a)!=1: continue
  baseline=float(a.mean_methylation_percent.iloc[0])
  for _,r in ss.iterrows(): rows.append(dict(dataset=name,position=pos,base=r.base,mean_methylation_percent=r.mean_methylation_percent,difference_from_A_pp=r.mean_methylation_percent-baseline,n_sites=r.site_count))
pos=pd.DataFrame(rows);pos.to_csv(O/'full_flanking_position_effects.csv',index=False)
# load cancer site-level, summarize context adjusted differences per patient and compartment
p=SOURCE_DIR/'GSE70090_oxBS_site_level_annotated.csv.gz'
df=pd.read_csv(p,usecols=['patient','status','sixmer','primary_compartment','methylation_percent','coverage','chrom','position0'])
print('site records',df.shape)
# patient, compartment baseline tumor-normal; all 256 motif x patient x compartment
base=df.groupby(['patient','primary_compartment','status']).agg(n=('methylation_percent','size'),mean=('methylation_percent','mean')).reset_index()
base=base.pivot(index=['patient','primary_compartment'],columns='status',values='mean')
base['compartment_delta']=base['Tumor']-base['Normal']
by=df.groupby(['patient','primary_compartment','sixmer','status']).agg(n=('methylation_percent','size'),mean=('methylation_percent','mean')).reset_index()
by=by.pivot(index=['patient','primary_compartment','sixmer'],columns='status',values=['n','mean'])
by.columns=[a+'_'+b for a,b in by.columns];by=by.reset_index()
by['delta_pp']=by['mean_Tumor']-by['mean_Normal']
by=by.merge(base['compartment_delta'].reset_index(),on=['patient','primary_compartment'])
by['residual_delta_pp']=by['delta_pp']-by['compartment_delta']
by.to_csv(O/'cancer_motif_compartment_patient_effects.csv',index=False)
# stable motifs based on 4 patients in a compartment min 30 each status
f=by[(by.n_Normal>=30)&(by.n_Tumor>=30)]
q=f.groupby(['primary_compartment','sixmer']).agg(n_patients=('patient','nunique'),mean_residual=('residual_delta_pp','mean'),sd_residual=('residual_delta_pp','std'),n_positive=('residual_delta_pp',lambda x:(x>0).sum()),n_negative=('residual_delta_pp',lambda x:(x<0).sum()),min_n_normal=('n_Normal','min'),min_n_tumor=('n_Tumor','min')).reset_index()
q=q[q.n_patients==4].sort_values(['primary_compartment','mean_residual'])
q.to_csv(O/'cancer_motif_compartment_four_patient_summary.csv',index=False)
print('retention examples',q[q.sixmer.isin(['CGCGCG','ATCGAT'])].to_string(index=False)[:4000])
# patient-level all compartment motif effect after subtracting patient global change
allbase=df.groupby(['patient','status']).methylation_percent.mean().unstack();allbase['delta']=allbase['Tumor']-allbase['Normal']
allmot=df.groupby(['patient','sixmer','status']).agg(n=('methylation_percent','size'),mean=('methylation_percent','mean')).reset_index().pivot(index=['patient','sixmer'],columns='status',values=['n','mean'])
allmot.columns=[a+'_'+b for a,b in allmot.columns];allmot=allmot.reset_index()
allmot=allmot.merge(allbase['delta'].rename('patient_delta'),on='patient');allmot['motif_delta']=allmot.mean_Tumor-allmot.mean_Normal;allmot['residual_delta']=allmot.motif_delta-allmot.patient_delta
allmot.to_csv(O/'cancer_motif_four_patient_global_adjusted.csv',index=False)
print(allmot[allmot.sixmer.isin(['CGCGCG','ATCGAT'])].to_string(index=False))
