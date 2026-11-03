"""Bounded exact-certificate pilot for joint local H,C,N data."""
import argparse,json,sys,time,warnings
from pathlib import Path
from collections import Counter
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from local_dual_lp.system import build,letter_hist
from model import span,reduced_basis,gram,nullspace,lift_coordinate,weight
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha
from scipy.optimize import linprog,OptimizeWarning
import system as exact
assert Path(exact.__file__).resolve()==(ROOT.parents[1]/'split/system.py').resolve()

# Every HiGHS call in this isolated process uses only one worker thread.
old_linprog=exact.linprog
warnings.filterwarnings('ignore',category=OptimizeWarning,message='Unrecognized options detected')
def one_thread(*args,**kw):
 kw['options']=dict(kw.get('options',{}),threads=1);return old_linprog(*args,**kw)
exact.linprog=one_thread

def main():
 p=argparse.ArgumentParser();p.add_argument('--max-tasks',type=int,default=12);args=p.parse_args()
 verified=json.loads((HERE/'verification.json').read_text());assert verified['status']=='PASS'
 for name,value in verified['hashes'].items():assert sha(Path(name))==value
 profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];Htypes=[]
 for r in profiles:
  if (r['q'],r['t'])!=(3,2):continue
  H=[normalize_word(w,r['p'],8) for w in r['H']];hist=tuple(sorted(letter_hist(span(H),8).items()))
  if hist not in Htypes:Htypes.append(hist)
 assert len(Htypes)==1
 source=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];groups={}
 for c in source:
  if c['representative']['hull']!=3:continue
  C=c['representative']['generators'];rad=reduced_basis([lift_coordinate(v,C) for v in nullspace(gram(C,9),5)]);ws=span(rad)
  planes={tuple(sorted((a,b,a^b))) for a in ws if a for b in ws if b and a!=b}
  for plane in planes:
   wt=[weight(w,9) for w in plane];support=[(w&511)|(w>>9) for w in plane]
   if sum(w<8 for w in wt)>2:continue
   if any(wt[a]==8 and any(b!=a and support[b]&~support[a]==0 for b in range(3)) for a in range(3)):continue
   Nhist=tuple(sorted(Counter([0]+wt).items()));Chist=tuple(sorted(Counter(weight(w,9) for w in span(C)).items()));key=(Chist,Nhist)
   groups.setdefault(key,[]).append(dict(cls=c['class_index'],N=list(plane[:2])))
 records=[]
 for (Chist,Nhist),ids in sorted(groups.items(),key=lambda kv:-len(kv[1]))[:args.max_tasks]:
  start=time.monotonic();E,b,G,V,labels=build(17,3,6,8,dict(Htypes[0]),dict(Chist),dict(Nhist),light_budget=2,weight8_min=120)
  # Scaling changes neither the feasible set nor the exact certificate system.
  scaleE=np.maximum(np.abs(E.astype(float)).max(axis=1),np.abs(b.astype(float)));scaleE=np.maximum(scaleE,1)
  scaleG=np.maximum(np.abs(G.astype(float)).max(axis=1),1)
  status=one_thread(np.zeros(len(V)),A_eq=E.astype(float)/scaleE[:,None],b_eq=b.astype(float)/scaleE,A_ub=G.astype(float)/scaleG[:,None],b_ub=np.zeros(len(G)),bounds=[(0,None)]*len(V),method='highs').status
  cert,msg=exact.farkas_exact(E,b,G) if status==2 else (None,'primal_status_'+str(status))
  record=dict(C_spectrum=dict(Chist),N_spectrum=dict(Nhist),cases=len(ids),sample=ids[:3],primal_status=status,status='EXCLUDED_EXACT_CELL' if cert else 'NO_CERTIFICATE',message=str(msg),seconds=time.monotonic()-start)
  if cert:
   nu,lam=cert;assert exact.verify(E,b,G,nu,lam);record['certificate']=dict(nu=list(map(str,nu)),lam=list(map(str,lam)))
  records.append(record);print(json.dumps({k:v for k,v in record.items() if k!='certificate'}),flush=True)
  atomic_json(HERE/'pilot.json',dict(records=records,total_groups=len(groups),processed=len(records),scope='Joint q3,t2 local data pilot; only exact certificates exclude their complete specified cells'))
if __name__=='__main__':main()
