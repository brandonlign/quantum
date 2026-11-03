"""Residual q1/t0, q2/t2 (hull1) and q2/t0 (hull3) exact fixed-H,N slices.

One TASK = (suffix class, (q,t) branch, N-orbit representative, distinguished
H class, radical-identification ORBIT representative, first-label chunk).
LIMIT/ERROR never exclude. Any FOUND is physically re-verified from scratch.
"""
import json,subprocess,sys,tempfile,time
from collections import Counter
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,lift_coordinate,reduced_basis,span,symp,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import payload
from kernel_context_cache.common import kernel_contexts,identification_orbits
from fixed_kernel.symmetry import suffix_slice
from kernel_y_cache.solve import solve
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

Z_BINARY=ROOT/'kernel_z_compact/z_enum_compact'
Z_JOINT=ROOT/'kernel_z_joint/z_enum_joint'
PLAN=ROOT/'kernel_case_plan/plan.json'
CHUNK=64  # first-label representatives per task chunk; part of the semantics
BRANCH_HULL={(1,0):1,(2,2):1,(2,0):3}

SEMANTIC=['rankone_residual/common.py','rankone_residual/search.py','rankone_residual/PROOF.md','kernel_case_plan/provenance.json',
 'kernel_context_cache/common.py','kernel_context_cache/PROOF.md','kernel_context_cache/verification.json',
 'kernel_z_compact/z_enum_compact.c','kernel_z_compact/z_enum_compact','kernel_z_compact/PROOF.md','kernel_z_compact/verification.json',
 'kernel_z_joint/z_enum_joint.c','kernel_z_joint/z_enum_joint','kernel_z_joint/build.py','kernel_z_joint/PROOF.md','kernel_z_joint/verification.json','kernel_z_joint/threshold_verification.json',
 'kernel_y_parity/solve.py','kernel_y_parity/PROOF.md','kernel_y_parity/verification.json',
 'kernel_y_cache/solve.py','kernel_y_cache/build.py','kernel_y_cache/PROOF.md','kernel_y_cache/verification.json',
 'fixed_kernel/common.py','fixed_kernel/symmetry.py','fixed_kernel/PROOF.md','fixed_kernel/symmetry_verification.json','fixed_kernel/joint_verification.json','fixed_kernel/generic_verification.json',
 'kernel_case_plan/plan.json','kernel_case_plan/PROOF.md','kernel_case_plan/prepare.py','kernel_capacity/profiles.json',
 'opus_parallel/suffix9_dimension5/suffix9_dim5.json','model.py','fixed_suffix_search.py','two_site.py','mixed_frames/frames.py',
 'm8r1/suffix_symmetries.py','gauge_rankone_h5/common.py','gauge_rankone_h5/symmetry.py','mixed_frames/stage3/common.py','r1_column/Q4_KERNEL.md']

def fingerprints():return {p:sha(ROOT/p) for p in SEMANTIC}

def validate(*,require_controls):
 plan=json.loads(PLAN.read_text());assert plan['status']=='COMPLETE_PLAN_RELATIVE_TO_CLASSIFICATIONS'
 for path,digest in plan['hashes'].items():assert sha(ROOT/path)==digest,path
 prov=json.loads((ROOT/'kernel_case_plan/provenance.json').read_text());assert prov['status']=='PASS_TRANSITIVE_SAVED_PROVENANCE'
 for path,digest in prov['hashes'].items():assert sha(ROOT/path)==digest,('provenance',path)
 assert prov['hashes'].get('kernel_case_plan/plan.json',sha(PLAN))==sha(PLAN)
 proofs=['kernel_z_compact/verification.json','kernel_z_joint/verification.json','kernel_z_joint/threshold_verification.json','kernel_context_cache/verification.json','kernel_y_parity/verification.json','kernel_y_cache/verification.json','fixed_kernel/symmetry_verification.json','fixed_kernel/joint_verification.json']
 if require_controls:proofs.append('rankone_residual/controls.json')
 for f in proofs:
  r=json.loads((ROOT/f).read_text());assert r['status']=='PASS',f
  for path,digest in r['hashes'].items():assert sha(ROOT/path)==digest,(f,path)
 assert Counter((p['q'],p['t']) for p in plan['prefix'])[(1,0)]==4
 return plan

def prefixes(plan,q,t):return [p for p in plan['prefix'] if (p['q'],p['t'])==(q,t)]

@lru_cache(maxsize=8)
def target(C):
 model=ExtensionModel(17,3,6,8,[255<<8],list(C));assert model.C_distance>=6;return model

@lru_cache(maxsize=64)
def slice_setup(C,N,profile_H,profile_p,q,t):
 """Contexts, identification orbits; cached by (C,N,H). Pure function of inputs."""
 model=target(C);N=list(N);data=quotient_input(model,N)
 assert model.hull==BRANCH_HULL[(q,t)] and len(reduced_basis(N))==q-1
 assert len(data['Gamma'])==14-q and len(reduced_basis(data['Gamma']))==14-2*q+t
 H=[normalize_word(h,profile_p,8) for h in profile_H]
 assert len(reduced_basis(H))==q+1 and min(weight(w,8) for w in span(H) if w)>=6
 H,order,contexts=kernel_contexts(H,8,data['Gamma'],with_group=True)
 assert all((c['g'],c['t'])==(q-t,t) for c in contexts)
 io=identification_orbits(H,8,contexts)
 members=sorted(i for o in io['orbits'] for i in o['members']);assert members==list(range(len(contexts)))
 return model,data,H,order,contexts,io

def z_search(model,data,c,reps,seconds,*,group=True,fixed=None,max_nodes=None,backend='compact'):
 assert backend in ('compact','joint')
 cmd=[str(Z_BINARY if backend=='compact' else Z_JOINT),'--gamma','--avoid-radical','--first-orbits','--kernel-data']
 if backend=='joint':cmd+=['--mitm4','65536']
 if seconds is not None:cmd+=['--seconds',str(seconds)]
 if max_nodes is not None:cmd+=['--max-nodes',str(max_nodes)]
 if fixed is not None:cmd+=['--prefix',','.join(map(str,fixed))]
 with tempfile.TemporaryDirectory() as directory:
  if group:
   g=Path(directory)/'prefix.txt';g.write_text(f'{len(c["group"])} {model.m-1}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');cmd+=['--prefix-group',str(g)]
  proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
  out,err=proc.communicate(payload(model,data,reps,c))
 assert proc.returncode==0,err
 stats=json.loads(err);leaves=[tuple(map(int,l.split()[1:])) for l in out.splitlines()]
 assert len(leaves)==stats['leaves']==len(set(leaves))
 return stats,leaves

def physical_witness(model,columns,q,d=6):
 """Independent of syndrome costs: explicit S,D ranks, isotropy, full 2^20 dual span.
 Exact N gives rank M = 14-q (NOT model.b)."""
 rank=len(reduced_basis(columns));crit=bool(model.criterion(columns))
 S,D=model.construct(columns);S=reduced_basis(S);D=reduced_basis(D)
 Swords=set(span(S))
 iso=all(symp(u,v,17)==0 for u in S for v in S)
 commute=all(symp(u,v,17)==0 for u in S for v in D)
 dmin=min(weight(w,17) for w in span(D) if w not in Swords)
 smin=min(weight(w,17) for w in Swords if w)
 ok=crit and rank==14-q and len(S)==14 and len(D)==20 and iso and commute and dmin>=d
 return dict(S=S,D=D,M_rank=rank,expected_M_rank=14-q,criterion=crit,S_rank=len(S),D_rank=len(D),isotropic=iso,S_in_Dperp=commute,
             logical_distance=dmin,stabilizer_min_weight=smin,target_distance=d,is_17_3_6=ok and d==6,valid_at_d=ok)

def run_task(row,branch,oi,prefix,ci,chunk,seconds,hashes,*,y_seconds=None,candidate_dir=None,backend='compact'):
 assert backend in ('compact','joint')
 q,t=branch['q'],branch['t'];N=branch['N_orbits'][oi]['N']
 model,data,H,order,contexts,io=slice_setup(tuple(row['C']),tuple(N),tuple(prefix['H']),prefix['p'],q,t)
 assert ci in [o['representative'] for o in io['orbits']]
 c=contexts[ci];ss=suffix_slice(model,N,data,c);reps=ss['first_representatives']
 nchunks=max(1,-(-len(reps)//CHUNK));assert 0<=chunk<nchunks;mine=reps[chunk*CHUNK:(chunk+1)*CHUNK]
 start=time.monotonic();stats,leaves=(z_search(model,data,c,mine,seconds,backend=backend) if mine else (dict(status='COMPLETE',nodes=0,leaves=0,reason=None),[]))
 evidence=[];tally=Counter();found=None
 for xs in leaves:
  a=solve(model,H,[lift_coordinate(x,data['units']) for x in xs],seconds=y_seconds)
  tally[a['status']+':'+str(a['reason'])]+=1;evidence.append(dict(X=list(xs),result=a))
  if a['status']=='FOUND':found=a;break
 complete=stats['status']=='COMPLETE' and stats.get('reason') is None and found is None and all(e['result']['status']=='COMPLETE' and e['result']['fixed_z_excluded'] for e in evidence) and len(evidence)==len(leaves)
 r=dict(cls=row['cls'],hull=row['hull'],q=q,t=t,orbit_index=oi,N=N,N_members=branch['N_orbits'][oi]['members'],profile_index=prefix['profile_index'],
  context_index=ci,identification_orbit=next(o['members'] for o in io['orbits'] if o['representative']==ci),identifications=len(contexts),
  chunk=chunk,chunks=nchunks,chunk_size=CHUNK,first=mine,first_representatives_total=len(reps),orbit_digest=digest(ss),
  H=H,site_order=order,backend=backend,z=stats,tally=dict(tally),leaves=evidence,
  status='FOUND' if found else ('EXCLUDED_SLICE_CHUNK' if complete else 'LIMIT'),seconds=time.monotonic()-start,hashes=hashes,
  scope='One chunk of one exact (C,N,H,identification-orbit) slice in residual exclusive rank-one column; not a target exclusion')
 if found:
  # Persist the raw candidate BEFORE any expensive or asserting check.
  cdir=Path(candidate_dir) if candidate_dir is not None else HERE/'candidates';cdir.mkdir(exist_ok=True)
  raw=cdir/('CANDIDATE_'+task_name(row['cls'],q,t,oi,prefix['profile_index'],ci,chunk))
  atomic_json(raw,dict(record=r,model_C=row['C'],R=[255<<8],columns=found['solutions'][0]['model_columns'],hashes=hashes))
  try:r['physical']=physical_witness(model,found['solutions'][0]['model_columns'],q)
  except Exception as e:r['physical_error']=repr(e)
  r['candidate_file']=str(raw.relative_to(ROOT));atomic_json(raw,dict(record=r,model_C=row['C'],R=[255<<8],columns=found['solutions'][0]['model_columns'],hashes=hashes))
 return r,ss

def digest(ss):
 import hashlib
 key=json.dumps(dict(eligible=ss['eligible'],orbits=ss['orbits'],reps=ss['first_representatives']),sort_keys=True)
 return hashlib.sha256(key.encode()).hexdigest()

def task_name(cls,q,t,oi,profile,ci,chunk):return f'{cls}_q{q}t{t}_{oi}_{profile}_{ci}_{chunk}.json'
