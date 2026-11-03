"""LIMIT-stream prefix check + repeated paired timing on the complete q2t0 root."""
import json,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
from gauge_rankone_h5.common import normalize_word
rows=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];out=[]
for hull,q,ri,secs,reps_n in [(1,1,0,'2',1),(1,2,8,'2',1),(3,2,4,'4',3)]:
 cls,C=next((i,r['representative']['generators']) for i,r in enumerate(rows) if r['representative']['hull']==hull)
 model=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))];data=quotient_input(model,[] if q==1 else rad[:q-1]);row=profiles[ri]
 H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
 reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
 with tempfile.TemporaryDirectory() as dd:
  g=Path(dd)/'g.txt';g.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
  base=['--gamma','--avoid-radical','--first-orbits','--kernel-data','--seconds',secs,'--prefix-group',str(g)]
  for k in range(reps_n):
   res={}
   for name,b,extra in [('frozen',ROOT/'fixed_kernel/z_enum_fixed',[]),('mitm4096',HERE/'z_enum_mitm4',['--mitm4','4096'])]:
    p=subprocess.run([str(b)]+base+extra,input=payload(model,data,reps,c),capture_output=True,text=True,check=True);res[name]=(p.stdout.splitlines(),json.loads(p.stderr))
   f,m=res['frozen'],res['mitm4096'];short,long=(f[0],m[0]) if len(f[0])<=len(m[0]) else (m[0],f[0])
   r=dict(profile=ri,repeat=k,frozen_status=f[1]['status'],mitm_status=m[1]['status'],frozen_seconds=f[1]['seconds'],mitm_seconds=m[1]['seconds'],frozen_leaves=len(f[0]),mitm_leaves=len(m[0]),shorter_is_prefix=long[:len(short)]==short)
   out.append(r);print(json.dumps(r),flush=True)
(HERE/'gate2_prefix.json').write_text(json.dumps(out,indent=1)+'\n')
