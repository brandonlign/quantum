"""Exercise isolated process dispatch on two already-complete target tasks."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from radical_projection.pool_run import task,pool_hashes
from radical_projection.run import fingerprints
from mixed_frames.stage3.common import atomic_json,sha

def main():
 plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text());row=plan['records'][0];indices=[4,6];expected=fingerprints();provenance=pool_hashes();previous={}
 for i in indices:
  r=json.loads((HERE/'runs'/f'{row["cls"]}_{i}.json').read_text());assert r['status']=='EXCLUDED_FIXED_N' and r['hashes']==expected;previous[i]=r
 with ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(task,[(row,i,expected,provenance) for i in indices]))
 for i,r in zip(indices,results):
  assert r['status']=='EXCLUDED_FIXED_N' and r['z']['status']=='COMPLETE'
  for key in ('nodes','leaves','gamma_pruned','projection_pruned'):assert r['z'][key]==previous[i]['z'][key]
 atomic_json(HERE/'pool_verification.json',dict(status='PASS',physical_target_tasks=2,matching_complete_tree_counts=True,pool_hashes=provenance,search_hashes=expected,verification_sha256=sha(Path(__file__))));print('PASS: two process workers reproduced complete target tree counts',flush=True)
if __name__=='__main__':main()
