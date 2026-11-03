"""Paired complete fixed-first roots (hull1 C0, four q1 forms, first 115), no brute check, 3 repeats."""
import hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import fixed_kernel_mitm4.gate1 as G
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,identification_orbits
from gauge_rankone_h5.common import normalize_word
classes=json.loads((G.ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((G.ROOT/'kernel_capacity/profiles.json').read_text())['records']
C=classes[0]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[]);out=[]
for ri in range(4):
 row=profiles[ri];H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
 reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
 g=HERE/'group_tmp.txt';g.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
 for k in range(3):
  o0,s0=G.run(G.FROZEN,model,data,c,['--prefix','115','--prefix-group',str(g)],reps);o1,s1=G.run(G.NEW,model,data,c,['--prefix','115','--prefix-group',str(g),'--mitm4','4096'],reps)
  assert s0['status']==s1['status']=='COMPLETE' and o0==o1
  r=dict(profile=ri,repeat=k,group=len(c['group']),frozen_nodes=s0['nodes'],mitm_nodes=s1['nodes'],frozen_seconds=s0['seconds'],mitm_seconds=s1['seconds'],tested=s1['mitm_tested'],pruned=s1['mitm_pruned'],skipped=s1['mitm_skipped'],leaves=s1['leaves'])
  out.append(r);print(json.dumps(r),flush=True)
g.unlink()
(HERE/'roots_timing.json').write_text(json.dumps(out,indent=1)+'\n')
