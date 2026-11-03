"""Paired complete target roots; exact streams and counters must match."""
import json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,identification_orbits
from gauge_rankone_h5.common import normalize_word
from kernel_z_compact.verify import pair
from mixed_frames.stage3.common import atomic_json,sha

def main():
 rows=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
 cls,C=next((i,r['representative']['generators']) for i,r in enumerate(rows) if r['representative']['hull']==1)
 model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[])
 profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];records=[]
 # Borrow one pool's cores only during these short sequential measurements.
 pool=int((ROOT/'radical_projection/pool.pid').read_text());command=subprocess.check_output(['ps','-p',str(pool),'-o','command='],text=True)
 assert '-m radical_projection.pool_run' in command and os.getpgid(pool)==pool
 os.killpg(pool,signal.SIGSTOP)
 try:
  for ri,row in enumerate(profiles):
   if (row['q'],row['t'])!=(1,0):continue
   H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True)
   c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
   for repeat in range(2):
    r=dict(cls=cls,profile_index=ri,repeat=repeat,group=len(c['group']),first=115,**pair(model,data,c,group=True,fixed=[115]))
    assert r['status']=='COMPLETE';records.append(r);print(json.dumps(r),flush=True)
 finally:os.killpg(pool,signal.SIGCONT)
 report=dict(status='PASS',scope='Four q1 H forms, one complete first-image root each, two paired repeats; no whole-case exclusion',records=records,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'verify.py',HERE/'z_enum_compact.c',HERE/'z_enum_compact',ROOT/'fixed_kernel/z_enum_fixed']})
 atomic_json(HERE/'benchmark.json',report)

if __name__=='__main__':main()
