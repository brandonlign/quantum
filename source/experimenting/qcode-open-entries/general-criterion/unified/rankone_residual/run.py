"""Bounded residual runner. No unlimited sweep without --all-classes and explicit
--max-tasks or primary approval. LIMIT records remain unresolved."""
import argparse,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from mixed_frames.stage3.common import atomic_json
from rankone_residual.common import validate,fingerprints
from rankone_residual.search import expected_slices,run_slice

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--classes',type=int,nargs='+',required=True)
 ap.add_argument('--branch',type=str,choices=['1,0','2,2','2,0'])
 ap.add_argument('--profiles',type=int,nargs='+')
 ap.add_argument('--seconds',type=float,help='Z time limit per chunk; LIMIT never excludes')
 ap.add_argument('--max-slices',type=int)
 ap.add_argument('--rewrite',action='store_true')
 ap.add_argument('--backend',choices=['compact','joint-auto'],default='compact')
 ap.add_argument('--smoke-without-controls',action='store_true',help='ONLY for building controls; writes to smoke_runs/')
 a=ap.parse_args();assert __debug__
 plan=validate(require_controls=not a.smoke_without_controls);hashes=fingerprints()
 folder=HERE/('smoke_runs' if a.smoke_without_controls else 'runs');folder.mkdir(exist_ok=True);done=0;summary=[]
 branches=None if a.branch is None else {tuple(map(int,a.branch.split(',')))}
 rows={r['cls']:r for r in plan['records']}
 for cls in a.classes:
  row=rows[cls]
  for s in expected_slices(plan,row,branches=branches,profiles=a.profiles):
   backend='joint' if a.backend=='joint-auto' and s['prefix']['profile_index'] in (0,8) else 'compact'
   st=run_slice(plan,row,s,a.seconds,hashes,folder,rewrite=a.rewrite,backend=backend)
   summary.append(dict(cls=cls,q=s['q'],t=s['t'],orbit=s['oi'],profile=s['prefix']['profile_index'],context=s['ci'],statuses=st))
   done+=1
   if a.max_slices is not None and done>=a.max_slices:break
  if a.max_slices is not None and done>=a.max_slices:break
 atomic_json(HERE/'finished.json',dict(status='BOUNDED_BATCH_FINISHED',args=vars(a),slices=summary,hashes=hashes))
if __name__=='__main__':main()
