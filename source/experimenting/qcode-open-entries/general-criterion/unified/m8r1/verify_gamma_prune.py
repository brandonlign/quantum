"""Check: on COMPLETE runs, gamma-pruned leaves == plain leaves filtered by z_condition (exact set)."""
import json,subprocess,sys,random
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from model import ExtensionModel
from run_benchmark import z_condition
ROOT=HERE.parents[2]
def run(binary,model,d,gamma,limit):
    inp=f'{model.m} {len(model.distances)} {d}\n'+' '.join(map(str,model.distances))+'\n'
    if gamma:inp+=' '.join(map(str,model.Gamma))+'\n'
    cmd=[str(HERE/binary),'--max-nodes',str(limit)]+(['--gamma'] if gamma else [])
    p=subprocess.run(cmd,input=inp,capture_output=True,text=True)
    st=json.loads(p.stderr.strip());L=[tuple(map(int,l.split()[1:])) for l in p.stdout.splitlines()]
    return st,L
cells=json.loads((ROOT/'release16/census/cells.json').read_text())
cases=[]
for n,m,key in [(16,6,'10,10'),(16,7,'9,8'),(16,8,'8,6'),(17,7,'9,8'),(17,8,'8,6')]:
    for i in range(len(cells[key])):cases.append((n,m,key,i))
random.seed(1);random.shuffle(cases)
done=agree=0;tot_plain=tot_g=0;out=[]
for n,m,key,i in cases:
    k=n-11 if n==16 else n-11
    try:model=ExtensionModel(n,k,5,m,[((1<<m)-1)<<m],cells[key][i])
    except AssertionError:continue
    s1,L1=run('z_enum8192',model,5,False,3_000_000)
    if s1['status']!='COMPLETE':continue
    s2,L2=run('z_enum_gamma',model,5,True,3_000_000)
    assert s2['status']=='COMPLETE'
    F=sorted(x for x in L1 if z_condition(x,model.Gamma,model.b))
    ok=F==sorted(L2);done+=1;agree+=ok
    out.append(dict(n=n,m=m,cell=key,cls=i,plain_nodes=s1['nodes'],gamma_nodes=s2['nodes'],plain_leaves=len(L1),filtered=len(F),gamma_leaves=len(L2),agree=ok))
    if not ok:print('DISAGREE',out[-1]);break
    if done>=60:break
print(json.dumps(dict(status='PASS' if agree==done and done else 'FAIL',instances=done,with_surviving_leaves=sum(o['filtered']>0 for o in out),nonempty_plain=sum(o['plain_leaves']>0 for o in out))))
json.dump(dict(instances=out),open(HERE/'gamma_prune_verification.json','w'),indent=1)
