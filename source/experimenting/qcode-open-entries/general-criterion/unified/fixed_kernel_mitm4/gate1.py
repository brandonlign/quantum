"""Gate 1: byte-identical COMPLETE LEAF streams, frozen native vs --mitm4.
Frozen and previous namespaces are only read/executed. Writes gate1.json."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace
from fixed_suffix_search import quotient_input
from gauge_rankone_h5.common import normalize_word
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
import fixed_kernel_lookahead.gate1 as prev
FROZEN=ROOT/'fixed_kernel/z_enum_fixed';NEW=HERE/'z_enum_mitm4'
BASE=['--gamma','--avoid-radical','--first-orbits','--kernel-data']
CAP=os.environ.get('CAP','4096');CHECK=os.environ.get('CHECK','2000000')
LEAF_TIME=25

def run(binary,model,data,c,extra,reps=None,timeout=LEAF_TIME):
 if reps is None:reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
 p=subprocess.run([str(binary)]+BASE+extra,input=payload(model,data,reps,c),capture_output=True,text=True,check=True,timeout=timeout)
 return p.stdout,json.loads(p.stderr)

def compare(name,model,data,c,extra,rows,reps=None,expect=None,need_complete=True,group='both'):
 modes=[]
 if group=='both':modes=[('g',True),('',False)]
 elif group=='only':modes=[('g',True)]
 else:modes=[('',False)]
 for mode,grouped in modes:
  real=list(extra)
  if grouped:
   g=HERE/'group_tmp.txt';g.write_text(f'{len(c["group"])} {model.m-1}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');real+=['--prefix-group',str(g)]
  out0,s0=run(FROZEN,model,data,c,real,reps)
  out1,s1=run(NEW,model,data,c,real+['--mitm4',CAP,'--mitm4-check',CHECK],reps)
  ok=out0==out1 and s0['status']==s1['status'] and (s0['status']=='COMPLETE' or not need_complete)
  if expect is not None:ok=ok and ('LEAF '+expect) in out1.splitlines()
  r=dict(name=name,g=c['g'],t=c['t'],group=grouped,status=s1['status'],leaves=s1['leaves'],frozen_nodes=s0['nodes'],nodes=s1['nodes'],frozen_seconds=s0['seconds'],seconds=s1['seconds'],
   **{k:s1[k] for k in ('mitm_tested','mitm_pruned','mitm_skipped','mitm_overflow','mitm_checked')},sha256=hashlib.sha256(out1.encode()).hexdigest(),identical=ok)
  rows.append(r);print(json.dumps({k:v for k,v in r.items() if k!='sha256'}),flush=True);assert ok,r

def large(rows):
 import random
 from model import coordinates,reduced_basis
 from fixed_kernel.common import standard
 from fixed_kernel.verify_large import paired_basis
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
 ref={r['hull']:r for r in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']}
 for hull,t in [(3,2),(5,0)]:
  # Same construction as fixed_kernel/verify_large.py (minus its 2^20 physical check).
  ci,C=next((i,r['representative']['generators']) for i,r in enumerate(classes) if r['representative']['hull']==hull)
  M6=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),5)];data=quotient_input(M6,rad[:2])
  prof=next(r for r in profiles if (r['q'],r['t'])==(3,t));H=[normalize_word(h,prof['p'],8) for h in prof['H']];H,_,cs=kernel_contexts(H,8,data['Gamma'])
  model=ExtensionModel(17,3,3,8,[255<<8],C);frame=standard(8);p=7;mask=127
  K=reduced_basis([coordinates(h,frame.R+frame.E)>>1 for h in H])
  sp=lambda a,b:(((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())%2
  A=nullspace([(k>>p)|((k&mask)<<p) for k in K],14);BA,g=paired_basis(A,sp);Gamma=data['Gamma'];b=len(Gamma)
  BQ,_=paired_basis([1<<j for j in range(b)],lambda a,c:(a&lift_coordinate(c,Gamma)).bit_count()%2)
  tau=[lift_coordinate(coordinates(1<<j,BQ),BA) for j in range(b)];rng=random.Random(1736+hull)
  for attempt in range(100):
   xs=[sum(sp(v,1<<i)<<j for j,v in enumerate(tau)) for i in range(14)];mc=frame.model_columns(model,[lift_coordinate(x,data['units']) for x in xs])
   if model.criterion(mc):break
   v=lift_coordinate(rng.randrange(1<<len(A)),A);tau=[w^(v if sp(w,v) else 0) for w in tau]
  assert mc==ref[hull]['known_columns'] and H==ref[hull]['H']
  ri=[lift_coordinate(lift_coordinate(r,tau),frame.E) for r in cs[0]['radical_coeff']]
  c=next(c for c in cs if all(a^b in (0,255<<8) for a,b in zip(c['radical_images'],ri)))
  for k in (7,6,5,4):
   compare(f'large17_h{hull}_prefix{k}',model,data,c,['--prefix',','.join(map(str,xs[:k]))],rows,expect=' '.join(map(str,xs[:7])),group='none')

def main():
 start=time.monotonic();rows=[]
 # small controls from previous namespace, with same compare contract (expect = leaf body).
 prev.compare=lambda name,model,data,c,extra,rs,need_complete=True,expect=None:compare(name,model,data,c,extra,rs,expect=(expect.split(' ',1)[1] if expect else None),need_complete=need_complete)
 prev.run=run
 prev.small(rows)
 large(rows)
 # Complete fixed-first target roots: hull1 class0, the four q1 H forms, first label 115.
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
 C=classes[0]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[])
 for ri in range(4):
  row=profiles[ri];H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
  reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
  first=115 if 115 in reps else reps[0]
  compare(f'target_h1_c0_prof{ri}_first{first}',model,data,c,['--prefix',str(first)],rows,reps=reps,group='only')
 report=dict(status='PASS',comparisons=len(rows),cap=int(CAP),brute_check_bound=int(CHECK),rows=rows,seconds=time.monotonic()-start,
  hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [FROZEN,ROOT/'fixed_kernel/z_enum_fixed.c',NEW,HERE/'z_enum_mitm4.c',HERE/'patch.py',Path(__file__)]})
 (HERE/'gate1.json').write_text(json.dumps(report,indent=1)+'\n');(HERE/'group_tmp.txt').unlink(missing_ok=True)
 print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
if __name__=='__main__':main()
