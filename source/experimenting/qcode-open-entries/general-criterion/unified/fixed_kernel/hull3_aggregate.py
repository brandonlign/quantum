"""Strict read-only aggregation of the q3,t2,hull3 extension branch."""
import json,sys
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from fixed_kernel.hull3_run import validate,fingerprints
from fixed_kernel.common import permute_word
from fixed_kernel.hull3_pool import pool_hashes
from model import reduced_basis,span,weight
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json

def main():
 plan=validate();hashes=fingerprints();classes=[];counts=Counter();totalH=Counter()
 for row in plan['records']:
  tally=Counter();errors=[]
  for oi,o in enumerate(row['orbits']):
   if o['screen'] is not None:tally['SCREENED']+=1;continue
   ht=Counter()
   for p in plan['prefix']:
    path=HERE/'hull3_runs'/f'{row["cls"]}_{oi}_{p["profile_index"]}.json'
    if not path.exists():ht['PENDING']+=1;continue
    try:
     r=json.loads(path.read_text());assert r['hashes']==hashes
     if 'pool_hashes' in r:assert r['pool_hashes']==pool_hashes()
     assert (r['cls'],r['orbit_index'],r['profile_index'])==(row['cls'],oi,p['profile_index']) and r['N']==o['N']
     expected=[normalize_word(h,p['p'],8) for h in p['H']]
     assert sorted(r['site_order'])==list(range(8))
     assert reduced_basis([permute_word(h,r['site_order'],8) for h in expected])==r['H']
     assert len(reduced_basis(r['H']))==4 and min(weight(w,8) for w in span(r['H']) if w)>=6
     if r['status']=='EXCLUDED_FIXED_H_N':
      assert r['z']['status']=='COMPLETE' and r['z']['reason'] is None and r['found'] is None
      assert r['z']['leaves']==len(r['leaves']) and sum(r['tally'].values())==len(r['leaves'])
      assert all(e['result']['status']=='COMPLETE' and e['result']['fixed_z_excluded'] and not e['result']['solutions'] for e in r['leaves'])
      assert r['context']['g']==1 and r['context']['t']==2
      assert len(r['eligible_first'])==len(set(r['eligible_first'])) and len(r['context']['radical_coeff'])==1
      data=r['orbit_data'];assert data['eligible']==r['eligible_first']
      members=[v for o in data['orbits'] for v in o['members']]
      assert len(members)==len(set(members)) and sorted(members)==data['eligible']
      assert data['first_representatives']==[min(o['members']) for o in data['orbits']]
      assert data['first_representatives']==[o['representative'] for o in data['orbits']]
      ht['EXCLUDED_FIXED_H_N']+=1
     else:ht[r['status']]+=1
    except Exception as e:ht['ERROR']+=1;errors.append(dict(orbit=oi,prefix=p['profile_index'],error=repr(e)))
   totalH.update(ht);tally['EXCLUDED_FIXED_N' if ht==Counter(EXCLUDED_FIXED_H_N=5) else 'OPEN_FIXED_N']+=1
  counts.update(tally);closed=set(tally)<=set(('SCREENED','EXCLUDED_FIXED_N'))
  classes.append(dict(cls=row['cls'],status='CONDITIONALLY_EXCLUDED_Q3_BRANCH' if closed else 'OPEN',counts=dict(tally),errors=errors))
 out=dict(status='CONDITIONALLY_COMPLETE' if all(c['status']=='CONDITIONALLY_EXCLUDED_Q3_BRANCH' for c in classes) else 'OPEN',counts=dict(counts),H_counts=dict(totalH),closed_classes=sum(c['status']=='CONDITIONALLY_EXCLUDED_Q3_BRANCH' for c in classes),classes=len(classes),records=classes,hashes=hashes,scope='ONLY q3,t2,hull3 exclusive rank-one column, conditional on complete C,N,H classifications and verified physical orbits; q2 hull3 and whole target still open')
 atomic_json(HERE/'hull3_aggregate.json',out);print(json.dumps({k:v for k,v in out.items() if k not in ('records','hashes')}),flush=True)
if __name__=='__main__':main()
