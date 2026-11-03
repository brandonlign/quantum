"""Exact full-context and outer-orbit equality on residual physical cases."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts as old,identification_orbits as old_orbits
from kernel_context_cache.common import kernel_contexts as new,identification_orbits as new_orbits,template,orbit_template
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def main():
 plan=json.loads((ROOT/'kernel_case_plan/plan.json').read_text());records=[]
 rows=[r for hull in [1,3] for r in [x for x in plan['records'] if x['hull']==hull][:2]]
 for row in rows:
  model=ExtensionModel(17,3,6,8,[65280],row['C'])
  for branch in row['branches']:
   q,t=branch['q'],branch['t'];N=branch['N_orbits'][0]['N'];data=quotient_input(model,N)
   for prefix in plan['prefix']:
    if (prefix['q'],prefix['t'])!=(q,t):continue
    H=[normalize_word(h,prefix['p'],8) for h in prefix['H']]
    start=time.monotonic();a=old(H,8,data['Gamma'],with_group=True);oa=old_orbits(a[0],8,a[2]);ot=time.monotonic()-start
    start=time.monotonic();b=new(H,8,data['Gamma'],with_group=True);ob=new_orbits(b[0],8,b[2]);nt=time.monotonic()-start
    assert a==b and oa==ob
    start=time.monotonic();c=new(H,8,data['Gamma'],with_group=True);oc=new_orbits(c[0],8,c[2]);warm=time.monotonic()-start
    assert c==a and oc==oa
    records.append(dict(cls=row['cls'],q=q,t=t,profile=prefix['profile_index'],contexts=len(b[2]),old_seconds=ot,new_seconds=nt,warm_seconds=warm));print(json.dumps(records[-1]),flush=True)
 report=dict(status='PASS',records=records,template_cache=template.cache_info()._asdict(),orbit_cache=orbit_template.cache_info()._asdict(),hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'PROOF.md',ROOT/'fixed_kernel/common.py',ROOT/'model.py',ROOT/'kernel_case_plan/plan.json']},scope='Exact returned data and orbit equality; prefix-only setup cache, no extension exclusion')
 atomic_json(HERE/'verification.json',report);print(json.dumps(dict(status='PASS',comparisons=len(records))),flush=True)

if __name__=='__main__':main()
