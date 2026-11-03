"""Bounded diagnostic: per-depth node counts of a COPY of z_enum_fixed (1 s each).
Mirrors fixed_kernel/pilot.py setup for one class per (q,t); never an exclusion."""
import json,subprocess,sys,tempfile
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
from gauge_rankone_h5.common import normalize_word

rows=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
out=[]
import os
EXTRA=os.environ.get('EXTRA','').split()
CASES=json.loads(os.environ.get('CASES','[[1,1,0,[0,3]],[1,2,2,[8]],[3,2,0,[4]]]'))
SECS=os.environ.get('SECS','1')
for hull,q,t,pis in CASES:
 cls,C=next((i,r['representative']['generators']) for i,r in enumerate(rows) if r['representative']['hull']==hull)
 model=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))]
 N=[] if q==1 else (rad[:1] if q==2 else rad[:2]);data=quotient_input(model,N)
 hist=Counter(data['distances'])
 for ri in pis:
  row=profiles[ri];H=[normalize_word(h,row['p'],8) for h in row['H']]
  H,order,contexts=kernel_contexts(H,8,data['Gamma'],with_group=True);ci=identification_orbits(H,8,contexts)['orbits'][0]['representative'];c=contexts[ci]
  reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
  with tempfile.TemporaryDirectory() as dd:
   p=Path(dd)/'g.txt';p.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
   r=subprocess.run([str(HERE/'z_enum_depth'),'--gamma','--avoid-radical','--first-orbits','--kernel-data','--seconds',SECS,'--prefix-group',str(p)]+EXTRA,input=payload(model,data,reps,c),capture_output=True,text=True)
  lines=r.stderr.strip().splitlines();stats=json.loads(lines[0]);depth=list(map(int,lines[1].split('LOOK')[0].split()[1:]))
  rec=dict(hull=hull,cls=cls,q=q,t=t,profile=ri,group=len(c['group']),labels=len(data['distances']),delta_hist=dict(sorted(hist.items())),first_reps=len(reps),status=stats['status'],nodes=stats['nodes'],leaves=stats['leaves'],prefix_pruned=stats['prefix_pruned'],depth_nodes=depth,look=lines[1].split('LOOK')[1].split(),extra=EXTRA)
  out.append(rec);print(json.dumps(rec),flush=True)
(HERE/os.environ.get('OUT','depth_probe.json')).write_text(json.dumps(dict(scope='1-second bounded diagnostics on a copied, instrumented enumerator; LIMIT is not exclusion',records=out),indent=1)+'\n')
