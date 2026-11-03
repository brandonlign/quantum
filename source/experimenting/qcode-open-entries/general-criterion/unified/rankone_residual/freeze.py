"""Freeze reviewed source and validation digests before the residual census."""
import argparse,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rankone_residual.common import validate,fingerprints
from mixed_frames.stage3.common import atomic_json,sha

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--backend-policy',choices=['compact','joint-auto'],default='compact');args=ap.parse_args();assert __debug__
 validate(require_controls=True);hashes=fingerprints()
 for name in ['controls.json','aggregation_verification.json','pool_verification.json']:
  path=HERE/name;r=json.loads(path.read_text());assert r['status']=='PASS',name
  for p,digest in r['hashes'].items():assert sha(ROOT/p)==digest,(name,p)
  hashes[str(path.relative_to(ROOT))]=sha(path)
 assert json.loads((HERE/'pool_verification.json').read_text())['backends']==['compact','joint']
 audits=['opus_parallel/independent_code_audit/COMPACT_AUDIT.md','opus_workers/reused_rankone/CACHE_AUDIT.md']
 if args.backend_policy=='joint-auto':audits.append('opus_parallel/kernel_audit/JOINT_AUDIT.md')
 for p in audits:assert (ROOT/p).exists(),p;hashes[p]=sha(ROOT/p)
 for p in [Path(__file__),HERE/'INTEGRATION.md',HERE/'pool.py',HERE/'aggregate.py',HERE/'verify_aggregate.py',HERE/'verify_pool.py',HERE/'controls.py',HERE/'run.py',ROOT/'scheduler/research_remaining.py',ROOT/'scheduler/rankone_after_888.py',ROOT/'scheduler/handoff_fast_frames.py']:
  hashes[str(p.relative_to(ROOT))]=sha(p)
 assert not list((HERE/'candidates').glob('*.json')),'Inspect saved candidates before starting'
 result=dict(status='READY_FOR_RESIDUAL_SWEEP',backend_policy=args.backend_policy,hashes=hashes,scope='Exact q1/t0 hull1, q2/t2 hull1 and q2/t0 hull3 only; both q3 branches and rank-two evidence remain separate')
 path=HERE/'READY.json'
 if path.exists():assert json.loads(path.read_text())==result,'A frozen version already exists; preserve its sources'
 else:atomic_json(path,result)
 print(json.dumps(dict(status=result['status'],backend_policy=args.backend_policy,pinned_files=len(hashes))),flush=True)

if __name__=='__main__':main()
