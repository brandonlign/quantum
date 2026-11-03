"""Bounded full-Pauli local LP pilots; exact certificates only."""
import json,sys,time,warnings
from pathlib import Path
from collections import Counter
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from pauli_orbit_lp.system import prepare,build
from pauli_orbit_lp.reduce import reduced,lift_certificate
from model import span
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha
from scipy.optimize import linprog,OptimizeWarning
import system as exact
warnings.filterwarnings('ignore',category=OptimizeWarning,message='Unrecognized options detected')
original=exact.linprog
def one_thread(*args,**kw):kw['options']=dict(kw.get('options',{}),threads=1,time_limit=30);return original(*args,**kw)
exact.linprog=one_thread

def main():
 verified=json.loads((HERE/'verification.json').read_text());assert verified['status']=='PASS'
 for name,value in verified['hashes'].items():assert sha(Path(name))==value
 profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];Hcases=[r for r in profiles if (r['q'],r['t'])==(3,2)]
 sample=next(r for r in json.loads((ROOT/'local_dual_lp/pilot.json').read_text())['records'] if r['N_spectrum'].get('6',0)==2);C_hist={int(w):v for w,v in sample['C_spectrum'].items()};N_hist={int(w):v for w,v in sample['N_spectrum'].items()};records=[]
 for row in Hcases:
  start=time.monotonic();H=[normalize_word(h,row['p'],8) for h in row['H']];data=prepare(H,8)
  print(json.dumps(dict(prepared=row['index'],group=data['physical_order'],A_orbits=len(data['AO']),B_orbits=len(data['BO']))),flush=True)
  if len(data['AO'])>200:print('SKIPPED_RESOURCE_BOUND',row['index'],flush=True);continue
  E,b,G,_=build(data,17,3,6,C_hist,N_hist,light_budget=2,weight8_min=120)
  er,br,gr,meta=reduced(E,b,G)
  print(json.dumps(dict(reduced_variables=er.shape[1],reduced_equations=len(br))),flush=True)
  scaleE=np.maximum(np.abs(er).max(axis=1),np.abs(br));scaleE=np.maximum(scaleE,1);scaleG=np.maximum(np.abs(gr).max(axis=1),1)
  result=one_thread(np.zeros(er.shape[1]),A_eq=er/scaleE[:,None],b_eq=br/scaleE,A_ub=gr/scaleG[:,None],b_ub=np.zeros(len(gr)),bounds=[(0,None)]*er.shape[1],method='highs-ipm')
  cert,message=exact.farkas_exact(er.astype(object),br.astype(object),gr.astype(object)) if result.status==2 else (None,result.message)
  record=dict(prefix_class=row['index'],H=H,C_hist=C_hist,N_hist=N_hist,A_orbits=len(data['AO']),B_orbits=len(data['BO']),variables=E.shape[1],equations=len(b),inequalities=len(G),primal_status=result.status,status='EXCLUDED_EXACT_CELL' if cert else 'NO_CERTIFICATE',message=str(message),seconds=time.monotonic()-start)
  if cert:
   nu,lam=lift_certificate(*cert,E,b,G,meta);assert exact.verify(E.astype(object),b.astype(object),G.astype(object),nu,lam);record['certificate']=dict(nu=list(map(str,nu)),lam=list(map(str,lam)))
  records.append(record);atomic_json(HERE/'pilot.json',dict(records=records,scope='One suffix/N spectrum across five valid q3,t2 prefix classes; only exact certificates exclude'))
  print(json.dumps({k:v for k,v in record.items() if k not in ('certificate','H')}),flush=True)
if __name__=='__main__':main()
