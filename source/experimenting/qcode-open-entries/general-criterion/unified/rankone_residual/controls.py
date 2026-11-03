"""Genuine [[17,3,3]] positives for every residual (q,t,hull) branch and every
distinguished H class, through the SAME Z backend, Y solver and witness check
used by production; plus identification/N coverage replay. Distance THREE only."""
import json,random,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,gram,lift_coordinate,nullspace,reduced_basis,span,split
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,standard
from fixed_kernel.verify_large import paired_basis
from fixed_kernel.symmetry import suffix_slice
from gauge_fix import supported_subspace
from rankone_residual.common import validate,fingerprints,slice_setup,z_search,solve,physical_witness,BRANCH_HULL
from mixed_frames.stage3.common import atomic_json,sha

def control(row,branch,prefix,rng):
 q,t=branch['q'],branch['t'];N=branch['N_orbits'][0]['N']
 M6,data,H,order,cs,io=slice_setup(tuple(row['C']),tuple(N),tuple(prefix['H']),prefix['p'],q,t)
 model=ExtensionModel(17,3,3,8,[255<<8],row['C']);frame=standard(8);p=7;mask=127
 K=reduced_basis([coordinates(h,frame.R+frame.E)>>1 for h in H]);assert len(K)==q
 sp=lambda a,b:(((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())%2
 A=nullspace([(k>>p)|((k&mask)<<p) for k in K],14)
 BA,g=paired_basis(A,sp);Gamma=data['Gamma'];b=len(Gamma)
 BQ,gq=paired_basis([1<<j for j in range(b)],lambda a,c:(a&lift_coordinate(c,Gamma)).bit_count()%2)
 assert g==gq==q-t and len(BA)==len(BQ)==b==14-q
 tau=[lift_coordinate(coordinates(1<<j,BQ),BA) for j in range(b)]
 assert gram([lift_coordinate(v,frame.E) for v in tau],8)==Gamma
 for attempt in range(400):
  xsmall=[sum(sp(v,1<<i)<<j for j,v in enumerate(tau)) for i in range(14)]
  columns=[lift_coordinate(x,data['units']) for x in xsmall];mc=frame.model_columns(model,columns)
  if model.criterion(mc):break
  v=lift_coordinate(rng.randrange(1<<len(A)),A);tau=[w^(v if sp(w,v) else 0) for w in tau]
 else:return dict(q=q,t=t,cls=row['cls'],profile=prefix['profile_index'],status='NO_D3_ISOMETRY_FOUND')
 rad=[lift_coordinate(lift_coordinate(r,tau),frame.E) for r in cs[0]['radical_coeff']]
 ci=next(i for i,c in enumerate(cs) if all(a^b in (0,255<<8) for a,b in zip(c['radical_images'],rad)))
 c=cs[ci];orbit=next(o for o in io['orbits'] if ci in o['members'])
 reps=[v for v,dd in enumerate(data['distances']) if dd>=model.d-1]
 pins={}
 for pin in (7,6,5):
  st,L=z_search(model,data,c,reps,20,group=False,fixed=xsmall[:pin]);assert st['status']=='COMPLETE' and tuple(xsmall[:7]) in L,(pin,st)
  sj,J=z_search(model,data,c,reps,20,group=False,fixed=xsmall[:pin],backend='joint')
  assert sj['status']=='COMPLETE' and J==L,(pin,sj)
  pins[pin]=len(L)
 ans=solve(model,H,columns[:7]);assert ans['status']=='FOUND'
 w=physical_witness(model,ans['solutions'][0]['model_columns'],q,d=3);assert w['valid_at_d'] and not w['is_17_3_6'],w
 w0=physical_witness(model,mc,q,d=3);assert w0['valid_at_d']
 S,D=model.construct(mc);actual=[split(x,8,9)[0] for x in supported_subspace(D,17,range(8))];assert reduced_basis(actual)==H
 return dict(q=q,t=t,hull=row['hull'],cls=row['cls'],profile=prefix['profile_index'],status='PASS',context=ci,identification_orbit=orbit['members'],pin_leaves=pins,
   M_rank=w['M_rank'],logical_distance=w['logical_distance'],physical_words_checked=1<<20,isometry_attempt=attempt+1)

def coverage(plan,row):
 out=[]
 for branch in row['branches']:
  q,t=branch['q'],branch['t'];assert row['hull']==BRANCH_HULL[(q,t)]
  mem=sorted(tuple(m) for o in branch['N_orbits'] for m in o['members'])
  if (q,t)==(2,0):assert len(mem)==7 and len(set(mem))==7
  else:assert len(mem)==1
  for prefix in [p for p in plan['prefix'] if (p['q'],p['t'])==(q,t)]:
   for oi,o in enumerate(branch['N_orbits']):
    model,data,H,order,cs,io=slice_setup(tuple(row['C']),tuple(o['N']),tuple(prefix['H']),prefix['p'],q,t)
    members=sorted(i for x in io['orbits'] for i in x['members']);assert members==list(range(len(cs)))==list(range((1,1,6)[q-t]))
    for x in io['orbits']:
     ss=suffix_slice(model,o['N'],data,cs[x['representative']])
     m=[v for y in ss['orbits'] for v in y['members']];assert sorted(m)==ss['eligible'] and len(m)==len(set(m))
    out.append(dict(q=q,t=t,orbit=oi,profile=prefix['profile_index'],identifications=len(cs),id_orbits=len(io['orbits'])))
 return out

def main():
 start=time.monotonic();plan=validate(require_controls=False);rng=random.Random(20261014);records=[];cov=[]
 rows=plan['records'];h1=[r for r in rows if r['hull']==1][:2];h3=[r for r in rows if r['hull']==3][:2]
 counts={(q,t):sum(1 for p in plan['prefix'] if (p['q'],p['t'])==(q,t)) for q,t in BRANCH_HULL};assert counts=={(1,0):4,(2,2):2,(2,0):6},counts
 for row in h1+h3:
  cov+=coverage(plan,row)
  for branch in row['branches']:
   for prefix in [p for p in plan['prefix'] if (p['q'],p['t'])==(branch['q'],branch['t'])]:
    r=control(row,branch,prefix,rng);records.append(r);print(json.dumps(r),flush=True)
 ok=[r for r in records if r['status']=='PASS'];branches={(r['q'],r['t']) for r in ok}
 status='PASS' if branches==set(BRANCH_HULL) and all(r['status']=='PASS' for r in records) else 'PARTIAL'
 report=dict(status=status,records=records,coverage=cov,seconds=time.monotonic()-start,
  scope='Distance-THREE positives and coverage replay only; not [[17,3,6]] constructions or exclusions',
  hashes={**fingerprints(),**{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'fixed_kernel/verify_large.py',ROOT/'gauge_fix.py']}})
 report['hashes'].pop('rankone_residual/controls.json',None)
 atomic_json(HERE/'controls.json',report);print(json.dumps(dict(status=status,controls=len(records),passed=len(ok),seconds=report['seconds'])),flush=True)
if __name__=='__main__':main()
