"""Adversarial proof-record checks on a genuinely complete target chunk."""
import copy,json,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rankone_residual.common import validate,fingerprints,run_task,slice_setup
from rankone_residual.search import expected_slices
from rankone_residual.aggregate import verify_chunk
from model import reduced_basis
from mixed_frames.stage3.common import atomic_json,sha

def main():
 plan=validate(require_controls=False);hashes=fingerprints();row=next(r for r in plan['records'] if r['cls']==5)
 s=next(expected_slices(plan,row,profiles=[4]));r,_=run_task(row,s['branch'],s['oi'],s['prefix'],s['ci'],0,None,hashes)
 assert r['status']=='EXCLUDED_SLICE_CHUNK'
 rejected=[]
 with tempfile.TemporaryDirectory() as directory:
  path=Path(directory)/s['names'][0];atomic_json(path,r);assert verify_chunk(path,row,s,0,hashes)=='EXCLUDED_SLICE_CHUNK'
  changes={
   'wrong_context':lambda x:x.update(context_index=x['context_index']+1),
   'wrong_first_chunk':lambda x:x.update(first=[]),
   'missing_Z_leaf':lambda x:x['z'].update(leaves=x['z']['leaves']+1),
   'Z_LIMIT_as_exclusion':lambda x:x['z'].update(status='LIMIT',reason='time_limit'),
   'wrong_identification_orbit':lambda x:x.update(identification_orbit=[]),
   'wrong_N_orbit':lambda x:x.update(N=[]),
   'stale_sources':lambda x:x.update(hashes={}),
   'unregistered_backend':lambda x:x.update(backend='unknown'),
  }
  model,data,*_=slice_setup(tuple(row['C']),tuple(r['N']),tuple(s['prefix']['H']),s['prefix']['p'],s['q'],s['t'])
  xs=[r['first'][0]]
  for i in range(len(data['Gamma'])):
   if len(reduced_basis(xs+[1<<i]))>len(xs):xs.append(1<<i)
   if len(xs)==7:break
  def fake_y(x,status):
   x['leaves']=[dict(X=xs,result=dict(status=status,complete=False,fixed_z_excluded=False,solutions=[],reason='search_limit',nodes=0))]
   x['z']['leaves']=1;x['tally']={status+':search_limit':1}
  changes['Y_LIMIT_as_exclusion']=lambda x:fake_y(x,'LIMIT')
  changes['Y_FOUND_as_exclusion']=lambda x:fake_y(x,'FOUND')
  for name,change in changes.items():
   bad=copy.deepcopy(r);change(bad);atomic_json(path,bad)
   try:verify_chunk(path,row,s,0,hashes)
   except (AssertionError,KeyError,ValueError):rejected.append(name)
   else:raise AssertionError('accepted corrupt record: '+name)
 report=dict(status='PASS',genuine_task=s['names'][0],complete_native_nodes=r['z']['nodes'],rejected=rejected,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'aggregate.py',HERE/'common.py',HERE/'search.py',HERE/'PROOF.md']})
 atomic_json(HERE/'aggregation_verification.json',report);print(json.dumps(report),flush=True)

if __name__=='__main__':main()
