import pandas as pd,numpy as np
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description="Sensitivity test for status-by-motif residuals")
parser.add_argument("--output-dir", default=".", help="Directory containing outputs from run_validation.py")
args=parser.parse_args()
P=Path(args.output_dir).resolve()
a=pd.read_csv(P/'site_summary_sufficient_statistics.csv')
out=[];mot=[]
for hold in [1,2,3,4]:
 tr=a[a.patient!=hold].copy();te=a[a.patient==hold].copy()
 b=tr.groupby('group').agg(n=('n','sum'),sy=('sy','sum'));b['mean']=b.sy/b.n
 tr['b']=tr.group.map(b['mean']);te['b']=te.group.map(b['mean'])
 tr['r']=tr.sy-tr.n*tr.b
 pooled=tr.groupby('sixmer').agg(n=('n','sum'),r=('r','sum'));pooled['effect']=pooled.r/(pooled.n+1000)
 tr['pooled']=tr.sixmer.map(pooled.effect)
 tr['r2']=tr.r-tr.n*tr.pooled
 ints=tr.groupby(['status','sixmer']).agg(n=('n','sum'),r=('r2','sum'));ints['effect']=ints.r/(ints.n+1500)
 te['pooled']=te.sixmer.map(pooled.effect)
 te['interaction']=[ints.loc[(s,m),'effect'] for s,m in zip(te.status,te.sixmer)]
 te['prediction']=(te.b+te.pooled+te.interaction).clip(0,1)
 for subset,t in [('all',te),('Normal',te[te.status=='Normal']),('Tumor',te[te.status=='Tumor'])]:
  sq=(t.sy2-2*t.prediction*t.sy+t.n*t.prediction**2).sum();n=t.n.sum()
  out.append({'patient':hold,'subset':subset,'RMSE_pp':100*np.sqrt(sq/n)})
 for status,t in te.groupby('status'):
  q=t.groupby('sixmer').apply(lambda z:pd.Series({'n':z.n.sum(),'pred':(z.n*z.prediction).sum()/z.n.sum()}),include_groups=False).reset_index()
  q['status']=status;q['patient']=hold;mot.append(q)
r=pd.DataFrame(out);r.to_csv(P/'heldout_status_interaction_metrics.csv',index=False)
m=pd.concat(mot).pivot(index=['patient','sixmer'],columns='status',values='pred');m['pred_delta_pp']=100*(m.Tumor-m.Normal)
d=pd.read_csv(P/'heldout_motif_tumor_delta.csv').set_index(['patient','sixmer'])
print(r.to_string(index=False))
print('motif tumor-normal delta RMSE',np.sqrt(np.mean((m.pred_delta_pp-d.observed_delta_pp)**2)))
