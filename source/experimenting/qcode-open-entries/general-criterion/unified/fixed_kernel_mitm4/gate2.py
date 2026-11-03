"""Gate 2: paired bounded target tests (frozen vs mitm4 caps) on the four
measured rows. SECS per run (default 2 s). LIMIT is never an exclusion."""
import json,os,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
from gauge_rankone_h5.common import normalize_word
SECS=os.environ.get('SECS','2')
rows=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];out=[]
for hull,q,t,ri in [(1,1,0,0),(1,1,0,3),(1,2,2,8),(3,2,0,4)]:
 cls,C=next((i,r['representative']['generators']) for i,r in enumerate(rows) if r['representative']['hull']==hull)
 model=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))]
 data=quotient_input(model,[] if q==1 else rad[:q-1]);row=profiles[ri]
 H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
 reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
 with tempfile.TemporaryDirectory() as dd:
  g=Path(dd)/'g.txt';g.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
  base=['--gamma','--avoid-radical','--first-orbits','--kernel-data','--seconds',SECS,'--prefix-group',str(g)]
  for name,binary,extra in [('frozen',ROOT/'fixed_kernel/z_enum_fixed',[]),('native_counting_copy',ROOT/'fixed_kernel_lookahead/z_enum_lookahead',[]),('mitm_cap4096',HERE/'z_enum_mitm4',['--mitm4','4096']),('mitm_cap65536',HERE/'z_enum_mitm4',['--mitm4','65536'])]:
   t0=time.monotonic();p=subprocess.run([str(binary)]+base+extra,input=payload(model,data,reps,c),capture_output=True,text=True,check=True)
   s=json.loads(p.stderr);r=dict(hull=hull,q=q,t=t,profile=ri,group=len(c['group']),first_reps=len(reps),run=name,status=s['status'],nodes=s['nodes'],leaves=s['leaves'],seconds=s['seconds'],wall=time.monotonic()-t0,
     depth_nodes=s.get('depth_nodes'),stdout_sha256=__import__('hashlib').sha256(p.stdout.encode()).hexdigest(),depth1_entered=(s.get('depth_nodes') or [None,None])[1],**{k:s.get(k) for k in ('mitm_tested','mitm_pruned','mitm_skipped','mitm_overflow')})
   out.append(r);print(json.dumps(r),flush=True)
(HERE/'gate2.json').write_text(json.dumps(dict(scope=f'{SECS}s bounded paired runs; LIMIT is not exclusion',records=out),indent=1)+'\n')
