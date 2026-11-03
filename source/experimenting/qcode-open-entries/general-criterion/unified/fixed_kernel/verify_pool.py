"""Two processes reproduce complete sequential native trees in isolated files."""
import json,sys,tempfile,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from fixed_kernel.hull3_run import validate,fingerprints
from fixed_kernel.hull3_pool import work
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();plan=validate();row=plan['records'][0];indices=[i for i,o in enumerate(row['orbits']) if o['screen'] is None][:2];hashes=fingerprints()
 with tempfile.TemporaryDirectory() as directory:
  tasks=[(row,i,plan['prefix'],hashes,None,directory,None) for i in indices]
  with ProcessPoolExecutor(max_workers=2) as pool:out=list(pool.map(work,tasks))
  for i in indices:
   with tempfile.TemporaryDirectory() as serial:
    work((row,i,plan['prefix'],hashes,None,serial,None))
    for p in plan['prefix']:
     name=f'{row["cls"]}_{i}_{p["profile_index"]}.json';a=json.loads((Path(directory)/name).read_text());b=json.loads((Path(serial)/name).read_text())
     assert a['status']==b['status']=='EXCLUDED_FIXED_H_N'
     assert a['leaves']==b['leaves'] and a['tally']==b['tally']
     assert all(a['z'][k]==b['z'][k] for k in ('status','nodes','leaves','gamma_pruned','xdom_pruned','prefix_pruned','projection_pruned','reason'))
 report=dict(status='PASS',workers=2,N_cases=len(indices),full_H_cases=10,independent_process_tree_counts_match=True,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'hull3_pool.py',HERE/'hull3_run.py',HERE/'pilot.py',HERE/'common.py',HERE/'solve.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed']})
 atomic_json(HERE/'hull3_pool_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
