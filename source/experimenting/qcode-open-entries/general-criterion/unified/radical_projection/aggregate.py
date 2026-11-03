"""Strict read-only aggregation of conditional hull-five coverage."""
import json,sys
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import reduced_basis,span,weight
from radical_projection.run import fingerprints,target
from radical_projection.pool_run import pool_hashes
from mixed_frames.stage3.common import atomic_json

def main():
 expected=fingerprints();plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text());screen=json.loads((ROOT/'rankone_recenter/column_screen.json').read_text());lookup={r['cls']:r for r in screen['records']};records=[];total=Counter()
 for row in plan['records']:
  tally=Counter();errors=[]
  for i,o in enumerate(row['orbits']):
   s=lookup[row['cls']]['orbits'][i];assert s['orbit_index']==i and s['covered_planes']==len(o['members'])
   N=o['representative'];assert len(reduced_basis(N))==2 and min(weight(w,9) for w in span(N) if w)>=6
   if s['reason'] is not None:tally['SCREENED']+=1;continue
   path=HERE/'runs'/f'{row["cls"]}_{i}.json'
   if not path.exists():tally['PENDING']+=1;continue
   try:
    r=json.loads(path.read_text());assert r['hashes']==expected
    if 'pool_hashes' in r:assert r['pool_hashes']==pool_hashes()
    assert r['cls']==row['cls'] and r['orbit_index']==i and r['N']==N
    if r['status']=='EXCLUDED_FIXED_N':
     assert r['z']['status']=='COMPLETE' and r['z']['reason'] is None
     assert r['z']['leaves']==len(r['leaves']) and r['found'] is None
     assert all(e['result']['status']=='EXCLUDED' for e in r['leaves'])
     assert sum(r['tally'].values())==len(r['leaves'])
     reps=r['orbit_data']['first_representatives'];assert len(reps)==len(set(reps))==r['orbit_data']['first_orbits']
     tally['EXCLUDED_FIXED_N']+=1
    else:tally[r['status']]+=1
   except Exception as e:tally['ERROR']+=1;errors.append(dict(orbit=i,error=str(e)))
  assert sum(tally.values())==len(row['orbits']);total.update(tally)
  closed=set(tally)<=set(('SCREENED','EXCLUDED_FIXED_N'))
  records.append(dict(cls=row['cls'],status='CONDITIONALLY_EXCLUDED_CLASS' if closed else 'OPEN',counts=dict(tally),errors=errors))
 result=dict(status='CONDITIONALLY_COMPLETE' if all(r['status']=='CONDITIONALLY_EXCLUDED_CLASS' for r in records) else 'OPEN',counts=dict(total),closed_classes=sum(r['status']=='CONDITIONALLY_EXCLUDED_CLASS' for r in records),classes=len(records),records=records,hashes=expected,scope='Hull-five exclusive rank-one column, conditional on suffix classification and recorded N/symmetry coverage; no whole-target exclusion')
 atomic_json(HERE/'aggregate.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('records','hashes')}))
if __name__=='__main__':main()
