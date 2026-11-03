"""Complete local H/R capacity profiles, conditional on checked small parents."""
import json,sys
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import reduced_basis,span,weight,symp,gram,dual,lift_coordinate
from classify_suffix9 import lift
from gauge_prefix_candidates import cert
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha
from mixed_frames.frames import xor_table

def profile(H,p,m=8):
 H=reduced_basis(H);assert p in span(H) and all(symp(p,h,m)==0 for h in H)
 seen=list(H);L=[]
 for w in dual([p],m):
  if len(reduced_basis(seen+[w]))>len(seen):seen.append(w);L.append(w)
 assert len(seen)==2*m-1
 hs=span(H);cost=[min(weight(w^h,m) for h in hs) for w in xor_table(L)]
 assert cost[0]==0 and all(w for w in cost[1:])
 return dict(histogram=dict(Counter(cost)),caps=[sum(c>=6-j for c in cost[1:]) for j in range(6)],quotient_bits=len(L))

def main():
 records=[]
 prefix=json.loads((ROOT/'gauge_prefix_candidates.json').read_text())['groups']
 for i,row in enumerate(prefix['2']['classes']):
  H=row['generators'];p=255<<8;records.append(dict(q=1,t=0,index=i,H=H,p=p,**profile(H,p)))
 parents=json.loads((ROOT/'opus_parallel/suffix9_dimension5/cells_d5_7_3.json').read_text())['7,3']['xz'];classes={};lifts=0
 # Independent unrestricted ordered-pair replay; local output GL is handled
 # solely by the physical certificate, not by a tail-orbit assumption.
 for g in parents:
  for a in range(8):
   for b in range(8):
    H=lift(g,a,b,7);ws=span(H)
    if min(weight(w,8) for w in ws if w)<6:continue
    full=[w for w in ws if weight(w,8)==8 and all(symp(w,h,8)==0 for h in H)]
    if not full:continue
    lifts+=1;key=cert(H,8)
    if key not in classes:classes[key]=dict(H=H,full=full)
 for i,row in enumerate(classes.values()):
  H=row['H'];t=len(reduced_basis(gram(H,8)));assert t in (0,2)
  # All distinguished words are kept; repeated capacity profiles are safe.
  for p in row['full']:
   normal=[normalize_word(h,p,8) for h in H];assert 255<<8 in span(normal)
   records.append(dict(q=2,t=t,index=i,H=H,p=p,**profile(H,p)))
 census=json.loads((ROOT/'gauge_rankone_h5/all_critical_dim4.json').read_text())['classes']
 for i,row in enumerate(census):
  H=row['generators'];t=row['gram_rank']
  for p in span(H):
   if weight(p,8)==8 and all(symp(p,h,8)==0 for h in H):records.append(dict(q=3,t=t,index=i,H=H,p=p,**profile(H,p)))
 report=dict(status='PASS_RELATIVE_TO_PARENT_CENSUS',records=records,dim3_classes=len(classes),dim3_ordered_lifts=lifts,groups={str(key):len([r for r in records if (r['q'],r['t'])==key]) for key in sorted({(r['q'],r['t']) for r in records})},hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'gauge_prefix_candidates.json',ROOT/'gauge_rankone_h5/all_critical_dim4.json',ROOT/'opus_parallel/suffix9_dimension5/cells_d5_7_3.json',ROOT/'model.py']},scope='Complete small H/R capacity candidates relative to saved parent censuses; no suffix or target exclusion')
 atomic_json(HERE/'profiles.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)
 for q,t in sorted({(r['q'],r['t']) for r in records}):print(q,t,sorted({tuple(r['caps']) for r in records if (r['q'],r['t'])==(q,t)}))
if __name__=='__main__':main()
