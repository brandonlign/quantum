"""Unlimited complete fixed-H sweep of hull2/hull6 888 cases by default."""
import argparse,json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from rank2_kernel.common import local_R,contexts,quotient_input,payload,identification_orbits
from rank2_kernel.symmetry import suffix_slice
from rank2_kernel.pilot import search
from rank2_kernel.verify_z import BINARY
from mixed_frames.stage3.common import atomic_json,sha

def validate():
 for name in ['prefix.json','suffix_plan.json','controls.json','y_verification.json','z_verification.json','symmetry_verification.json']:
  report=json.loads((HERE/name).read_text());assert report['status'].startswith('PASS'),name
  for path,digest in report['hashes'].items():assert sha(ROOT/path)==digest,path
 return json.loads((HERE/'suffix_plan.json').read_text()),json.loads((HERE/'prefix.json').read_text())['prefix']

def fingerprints():
 paths=[Path(__file__),HERE/'PROOF.md',HERE/'PREFIX.md',HERE/'common.py',HERE/'symmetry.py',HERE/'solve.py',HERE/'pilot.py',HERE/'verify_z.py',BINARY,BINARY.with_suffix('.c'),ROOT/'model.py',ROOT/'two_site.py',ROOT/'mixed_frames/frames.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'gauge_rankone_h5/symmetry.py',ROOT/'gauge_rankone_h5/common.py',ROOT/'frame_selection/certs.json',ROOT/'weight8_split/count_branches.json',ROOT/'opus_workers/rank2_backend/AUDIT.md']+[HERE/n for n in ['prefix.json','suffix_plan.json','controls.json','y_verification.json','z_verification.json','symmetry_verification.json']]
 return {str(p.relative_to(ROOT)):sha(p) for p in paths}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--classes',type=int,nargs='+');parser.add_argument('--q',type=int,nargs='+',choices=[0,1,3],default=[1,3]);parser.add_argument('--max-N',type=int);parser.add_argument('--seconds',type=float);args=parser.parse_args();assert __debug__
 plan,prefix=validate();hashes=fingerprints();folder=HERE/'runs';folder.mkdir(exist_ok=True);done=0
 for row in plan['records']:
  if row['q'] not in args.q or args.classes is not None and row['cls'] not in args.classes:continue
  model=None
  for oi,o in enumerate(row['orbits']):
   if o['screen']:continue
   if model is None:model=ExtensionModel(17,3,6,8,local_R(8),row['C']);assert model.hull==row['hull'] and model.C_distance>=6
   data=quotient_input(model,o['N'])
   for hi,p in enumerate(prefix[str(row['q'])]):
    H=p['H'];frame,cs=contexts(H,8,data['Gamma'],with_group=True);outer=identification_orbits(H,frame,cs)
    for orbit in outer['orbits']:
     ci=orbit['representative'];path=folder/f'{row["cls"]}_{oi}_{hi}_{ci}.json'
     if path.exists():
      old=json.loads(path.read_text());assert old['hashes']==hashes,path
      if old['status']=='EXCLUDED_IDENTIFICATION_SLICE':continue
      if old['status']=='FOUND':raise RuntimeError(f'Inspect target candidate {path}')
     c=cs[ci];ss=suffix_slice(model,o['N'],frame,data,c);start=time.monotonic();r=search(model,H,frame,data,c,ss,args.seconds)
     r.update(cls=row['cls'],q=row['q'],hull=row['hull'],N_orbit=oi,N=o['N'],H_index=hi,H=H,context_index=ci,context=c,identification_orbits=outer,orbit_data=ss,R=frame.R,E=frame.E,costs=frame.costs,Gamma=data['Gamma'],input_sha256=__import__('hashlib').sha256(payload(model,frame,data,c,ss['first_representatives']).encode()).hexdigest(),hashes=hashes,seconds=time.monotonic()-start,scope='One exact H,N and radical-identification orbit in 888; no whole-target exclusion')
     atomic_json(path,r);print(json.dumps({k:r[k] for k in ('cls','q','N_orbit','H_index','context_index','status','seconds','z')}),flush=True)
     if r['found']:raise RuntimeError(f'Pure target candidate saved in {path}')
   done+=1
   if args.max_N is not None and done>=args.max_N:return
 atomic_json(HERE/'finished.json',dict(status='SELECTED_BATCH_FINISHED' if args.classes is not None or set(args.q)!={1,3} else 'FINISHED_REQUIRES_STRICT_AGGREGATION',hashes=hashes))
if __name__=='__main__':main()
