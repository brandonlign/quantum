"""Replay every valid tail form and compare all full physical certificates."""
import json,sys,time
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import reduced_basis,span,weight,gram
from mixed_frames.frames import xor_table
from suffix_mass_check.common import canonical
from mixed_frames.stage3.common import atomic_json,sha
DATA=ROOT/'opus_parallel/suffix9_dimension5'

def main():
 start=time.monotonic();parent_check=json.loads((HERE/'verification8.json').read_text());assert parent_check['status']=='PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS'
 for name,value in parent_check['hashes'].items():assert sha(Path(name))==value
 parent=DATA/'cell_8_5.json';source=DATA/'suffix9_dim5.json';frozen={str(p):sha(p) for p in [parent,source]}
 parents=json.loads(parent.read_text())['xz'];classes=json.loads(source.read_text());expected={};multiplicity=[]
 for c in classes['classes']:
  g=c['representative']['generators'];assert len(reduced_basis(g))==5 and min(weight(w,9) for w in span(g) if w)>=6
  assert 5-len(reduced_basis(gram(g,9)))==c['representative']['hull']
  key=canonical(g,9);assert key not in expected;expected[key]=c['class_index'];multiplicity.append(len(c['extension_indices']))
 assert len(expected)==32865
 forms=[(0,0)]+[(a,0) for a in range(1,32)]+[(a,b) for a in range(1,32) for b in range(a+1,32) if b<a^b];assert len(forms)==187
 def zero_mask(a,b):return sum(int(not ((a&c).bit_count()%2 or (b&c).bit_count()%2))<<c for c in range(32))
 zeros=[zero_mask(a,b) for a,b in forms];allzeros=[zero_mask(a,b) for a in range(32) for b in range(32)];counts=[0]*len(expected);canonical_count=ordered_count=0;hulls=Counter()
 print('all 32865 source classes independently certified',flush=True)
 for pi,g in enumerate(parents):
  words=xor_table(g);assert len(set(words))==32;bad=sum(1<<c for c,w in enumerate(words) if weight(w,8)==5)
  ordered_count+=sum(not (bad&z) for z in allzeros)
  for (a,b),zero in zip(forms,zeros):
   if bad&zero:continue
   lifted=[(v&255)|((v>>8)<<9)|(((a>>i)&1)<<8)|(((b>>i)&1)<<17) for i,v in enumerate(g)]
   assert min(weight(w,9) for w in span(lifted) if w)>=6
   key=canonical(lifted,9);assert key in expected,'A valid lift is missing from the saved nine-site classification'
   counts[expected[key]]+=1;canonical_count+=1;hulls[5-len(reduced_basis(gram(lifted,9)))]+=1
  if pi%3000==0:print(json.dumps(dict(parents=pi,total_parents=len(parents),valid_lifts=canonical_count)),flush=True)
 assert counts==multiplicity and all(counts)
 assert canonical_count==classes['extension_representatives']==299175 and ordered_count==1791031
 assert dict(hulls)=={int(k):v for k,v in classes['hull_counts_reps'].items()}
 assert frozen=={name:sha(Path(name)) for name in frozen}
 report=dict(status='PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS',parent_verification_sha256=sha(HERE/'verification8.json'),parent_classes=len(parents),unrestricted_ordered_pairs=ordered_count,canonical_lifts=canonical_count,source_classes=len(expected),every_star_certificate_and_class_multiplicity_matches=True,hull_counts_reps=dict(hulls),class_multiplicities=counts,seconds=time.monotonic()-start,hashes={str(p):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'verification8.json',parent,source,ROOT/'model.py']},scope='Independent site-star encoding of all valid lifts; complete nine-site dimension-five distance-six coverage inherits the independently checked eight-site mass and certified release16 seven-site parent cells; no quantum extension exclusion')
 atomic_json(HERE/'verification9.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('class_multiplicities','hashes')}),flush=True)
if __name__=='__main__':main()
