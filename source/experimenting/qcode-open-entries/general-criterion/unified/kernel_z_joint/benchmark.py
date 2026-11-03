"""Bounded paired whole-root tests; LIMIT is never an exclusion."""
import hashlib,json,os,signal,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
from fixed_kernel.symmetry import suffix_slice
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def main():
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];records=[]
 pool=int((ROOT/'radical_projection/pool.pid').read_text());assert os.getpgid(pool)==pool
 assert '-m radical_projection.pool_run' in subprocess.check_output(['ps','-p',str(pool),'-o','command='],text=True)
 os.killpg(pool,signal.SIGSTOP)
 try:
  for hull,ri in [(1,0),(1,3),(1,8),(3,4)]:
   cls,C=next((i,r['representative']['generators']) for i,r in enumerate(classes) if r['representative']['hull']==hull);row=profiles[ri]
   model=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))];N=[] if row['q']==1 else rad[:row['q']-1];data=quotient_input(model,N)
   H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']];reps=suffix_slice(model,N,data,c)['first_representatives']
   with tempfile.TemporaryDirectory() as directory:
    g=Path(directory)/'group';g.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');text=payload(model,data,reps,c)
    base=[str(HERE/'z_enum_joint'),'--gamma','--avoid-radical','--first-orbits','--kernel-data','--prefix-group',str(g),'--seconds','2']
    for repeat in range(2):
     outputs=[]
     for mode,args in [('compact',[]),('joint',['--mitm4','65536'])]:
      start=time.monotonic();p=subprocess.run(base+args,input=text,capture_output=True,text=True,check=True,timeout=15);stats=json.loads(p.stderr);outputs.append(p.stdout)
      r=dict(cls=cls,profile_index=ri,repeat=repeat,mode=mode,first_reps=len(reps),group=len(c['group']),wall=time.monotonic()-start,stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest(),**stats)
      records.append(r);print(json.dumps(r),flush=True)
     assert outputs[0].startswith(outputs[1]) or outputs[1].startswith(outputs[0]),'different complete leaf ordering'
     if all(r['status']=='COMPLETE' for r in records[-2:]):assert outputs[0]==outputs[1]
 finally:os.killpg(pool,signal.SIGCONT)
 atomic_json(HERE/'benchmark.json',dict(status='PASS_BOUNDED_COMPARISON',scope='Four one-class one-context rows, two paired repeats at two CPU seconds; LIMIT is unresolved',records=records,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'z_enum_joint.c',HERE/'z_enum_joint',HERE/'verification.json']}))

if __name__=='__main__':main()
