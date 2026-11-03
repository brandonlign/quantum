"""Exact five-H sweep of the q3,t2,hull3 branch; limits never exclude."""
import argparse,json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,lift_coordinate,reduced_basis,span,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts
from fixed_kernel.pilot import search
from fixed_kernel.symmetry import suffix_slice
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def validate():
 plan=json.loads((HERE/'hull3_plan.json').read_text());assert plan['status']=='PASS_RELATIVE_TO_CLASSIFICATIONS'
 for path,digest in plan['hashes'].items():assert sha(ROOT/path)==digest,path
 for filename in ['generic_verification.json','large_verification.json','symmetry_verification.json','joint_verification.json']:
  proof=json.loads((HERE/filename).read_text());assert proof['status']=='PASS'
  for path,digest in proof['hashes'].items():assert sha(ROOT/path)==digest,path
 return plan

def fingerprints():
 paths=[Path(__file__),HERE/'pilot.py',HERE/'common.py',HERE/'solve.py',HERE/'symmetry.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',HERE/'PROOF.md',HERE/'hull3_plan.json',HERE/'generic_verification.json',HERE/'large_verification.json',HERE/'symmetry_verification.json',HERE/'joint_verification.json',ROOT/'fixed_suffix_search.py',ROOT/'two_site.py',ROOT/'model.py',ROOT/'mixed_frames/frames.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'gauge_rankone_h5/common.py',ROOT/'gauge_rankone_h5/symmetry.py',ROOT/'mixed_frames/stage3/common.py',ROOT/'opus_workers/fixed_kernel_final_audit/AUDIT.md']
 return {str(p.relative_to(ROOT)):sha(p) for p in paths}

def run_one(model,N,prefix,data,seconds,hashes,cls,orbit_index):
 H=[normalize_word(h,prefix['p'],8) for h in prefix['H']];assert min(weight(w,8) for w in span(H) if w)>=6
 H,order,contexts=kernel_contexts(H,8,data['Gamma'],with_group=True);assert len(contexts)==1
 c=contexts[0];assert c['g']==1 and c['t']==2
 orbits=suffix_slice(model,N,data,c);eligible=orbits['eligible'];reps=orbits['first_representatives']
 start=time.monotonic();r=search(model,H,data,c,reps,seconds)
 complete=r['z']['status']=='COMPLETE' and r['z']['reason'] is None and r['found'] is None and all(e['result']['status']=='COMPLETE' and e['result']['fixed_z_excluded'] for e in r['leaves'])
 assert r['z']['leaves']==len(r['leaves']) and sum(r['tally'].values())==len(r['leaves'])
 r.update(status='FOUND' if r['found'] else ('EXCLUDED_FIXED_H_N' if complete else 'LIMIT'),cls=cls,orbit_index=orbit_index,N=N,H=H,profile_index=prefix['profile_index'],site_order=order,context=c,eligible_first=eligible,orbit_data=orbits,seconds=time.monotonic()-start,hashes=hashes,scope='One full H and one exact N in q3,t2,hull3 exclusive rank-one column')
 if r['found']:
  candidate=r['found']['solutions'][0]['model_columns'];assert model.criterion(candidate);S,D=model.construct(candidate);assert min(weight(w,17) for w in span(D) if w)>=6
  assert len(reduced_basis(candidate))==11;r['physical_candidate']=dict(S=S,D=D)
 return r

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float);parser.add_argument('--max-tasks',type=int);parser.add_argument('--classes',type=int,nargs='+');args=parser.parse_args();assert __debug__
 plan=validate();hashes=fingerprints();folder=HERE/'hull3_runs';folder.mkdir(exist_ok=True);tasks=0
 for row in plan['records']:
  if args.classes is not None and row['cls'] not in args.classes:continue
  model=None
  for oi,orbit in enumerate(row['orbits']):
   if orbit['screen'] is not None:continue
   if model is None:model=ExtensionModel(17,3,6,8,[255<<8],row['C']);assert model.hull==3 and model.C_distance>=6
   N=orbit['N'];data=quotient_input(model,N)
   for prefix in plan['prefix']:
    path=folder/f'{row["cls"]}_{oi}_{prefix["profile_index"]}.json'
    if path.exists():
     prev=json.loads(path.read_text());assert prev['hashes']==hashes,path
     if prev['status']=='EXCLUDED_FIXED_H_N':continue
     if prev['status']=='FOUND':raise RuntimeError(f'Inspect target candidate {path}')
    r=run_one(model,N,prefix,data,args.seconds,hashes,row['cls'],oi);atomic_json(path,r);tasks+=1
    print(json.dumps({k:r[k] for k in ('cls','orbit_index','profile_index','status','seconds','z','tally')}),flush=True)
    if r['status']=='FOUND' or args.max_tasks is not None and tasks>=args.max_tasks:return
 atomic_json(HERE/'hull3_finished.json',dict(status='SELECTED_BATCH_FINISHED' if args.classes is not None else 'FINISHED_REQUIRES_STRICT_AGGREGATION',hashes=hashes))
if __name__=='__main__':main()
