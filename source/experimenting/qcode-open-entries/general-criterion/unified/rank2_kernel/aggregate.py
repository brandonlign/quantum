"""Strict read-only aggregation of the hull2/hull6 exact 888 kernel sweep."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,reduced_basis
from rank2_kernel.run import validate,fingerprints
from rank2_kernel.common import local_R,contexts,quotient_input,payload,identification_orbits
from rank2_kernel.symmetry import suffix_slice
from mixed_frames.stage3.common import atomic_json,sha

def normalized(v):return json.loads(json.dumps(v))

def main():
 plan,prefix=validate();hashes=fingerprints();counts=Counter();records=[];leaf_count=0;record_hashes={}
 for row in plan['records']:
  if row['q'] not in (1,3):continue
  tally=Counter();errors=[];model=None
  for oi,o in enumerate(row['orbits']):
   if o['screen'] is not None:tally['SCREENED_N']+=1;continue
   if model is None:model=ExtensionModel(17,3,6,8,local_R(8),row['C']);assert model.hull==row['hull'] and model.C_distance>=6
   data=quotient_input(model,o['N']);N_closed=True
   for hi,p in enumerate(prefix[str(row['q'])]):
    H=p['H'];frame,cs=contexts(H,8,data['Gamma'],with_group=True);outer=identification_orbits(H,frame,cs)
    for orbit in outer['orbits']:
     ci=orbit['representative'];path=HERE/'runs'/f'{row["cls"]}_{oi}_{hi}_{ci}.json'
     if not path.exists():tally['PENDING_SLICE']+=1;N_closed=False;continue
     try:
      r=json.loads(path.read_text());assert r['hashes']==hashes
      record_hashes[str(path.relative_to(ROOT))]=sha(path)
      assert (r['cls'],r['q'],r['hull'],r['N_orbit'],r['H_index'],r['context_index'])==(row['cls'],row['q'],row['hull'],oi,hi,ci)
      assert r['N']==o['N'] and r['H']==H and r['R']==frame.R and r['E']==frame.E and r['costs']==frame.costs and r['Gamma']==data['Gamma']
      assert r['context']==normalized(cs[ci]) and r['identification_orbits']==normalized(outer)
      ss=suffix_slice(model,o['N'],frame,data,cs[ci]);assert r['orbit_data']==normalized(ss)
      assert r['input_sha256']==hashlib.sha256(payload(model,frame,data,cs[ci],ss['first_representatives']).encode()).hexdigest()
      if r['status']=='EXCLUDED_IDENTIFICATION_SLICE':
       assert r['found'] is None and r['z']['status']=='COMPLETE' and r['z']['reason'] is None
       Xs=[tuple(e['X']) for e in r['leaves']];assert len(Xs)==len(set(Xs))==r['z']['leaves']
       assert all(len(x)==frame.p and len(reduced_basis(x))==frame.p and x[0] in ss['first_representatives'] for x in Xs)
       assert all(e['result']['status']=='COMPLETE' and e['result']['complete'] and e['result']['fixed_z_excluded'] and not e['result']['solutions'] for e in r['leaves'])
       tally['EXCLUDED_SLICE']+=1;leaf_count+=len(Xs)
      else:tally[r['status']]+=1;N_closed=False
     except Exception as e:tally['ERROR']+=1;N_closed=False;errors.append(dict(N_orbit=oi,H=hi,context=ci,error=repr(e)))
   tally['EXCLUDED_N' if N_closed else 'OPEN_N']+=1
  closed=tally['OPEN_N']==0
  records.append(dict(cls=row['cls'],q=row['q'],status='CONDITIONALLY_EXCLUDED_888_CLASS' if closed else 'OPEN',counts=dict(tally),errors=errors));counts.update(tally)
 result=dict(status='CONDITIONALLY_COMPLETE' if all(r['status']=='CONDITIONALLY_EXCLUDED_888_CLASS' for r in records) else 'OPEN',classes=len(records),closed_classes=sum(r['status']=='CONDITIONALLY_EXCLUDED_888_CLASS' for r in records),counts=dict(counts),Y_leaves_completely_rejected=leaf_count,records=records,hashes=hashes,record_hashes=record_hashes,aggregation_source_sha256=sha(Path(__file__)),scope='ONLY hull2/hull6 in exact rank-two 888 frame, conditional on complete prefix/suffix classifications; hull0 and rank-one columns remain separate')
 atomic_json(HERE/'aggregate.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('records','hashes','record_hashes')}),flush=True)
if __name__=='__main__':main()
