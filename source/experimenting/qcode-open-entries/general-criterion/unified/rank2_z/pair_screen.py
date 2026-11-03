"""Independent necessary unary/pair screen for nondegenerate rank-two extensions.

PASS is not FOUND; EXCLUDED applies only to one fixed C and fixed Z map.
"""
from pathlib import Path
from collections import Counter
import sys,json,time,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import reduced_basis
BASE=Path(__file__).resolve().parent


def linear_values(gens):
    out=[0]
    for g in gens:out += [v^g for v in out]
    return np.array(out,dtype=np.int64)


def cost_vector(m,x):
    p=m-2;mask=(1<<m)-1;out=[]
    for z in range(1<<p):
        physical_z=z|((z.bit_count()%2)<<(m-1))
        out.append(min(((x^rx)|(physical_z^rz)).bit_count() for rx in [0,mask] for rz in [0,mask]))
    return np.array(out,dtype=np.int64)


def screen(data,A,d=6,minimum=8):
    m=data.get('m',8);p=m-2;N=len(data['flags']);V=1<<p
    if len(A)!=p or len(reduced_basis(A))!=p or N!=1<<(2*p):raise ValueError('Require full-rank Lagrangian Z map in equal-dimensional nondegenerate spaces')
    form=data['form_rows'];delta=np.array(data['coset_distances']);lam=np.array(data['suffix_weights'])
    if any((form[a]&b).bit_count()%2 for a in A for b in A):raise ValueError('Z map is not isotropic')
    Az=linear_values(A);images=np.arange(N,dtype=np.int64)[:,None]^Az[None,:]
    def allowed(x):
        costs=cost_vector(m,x)
        return np.all((delta[images]>=d-costs[None,:]) & (lam[images]>=minimum-costs[None,:]),axis=1)
    # e=0 is separately represented by C and R; the zero prefix translate
    # should be omitted rather than tested against d or minimum.
    costs=cost_vector(m,0)
    if any(delta[int(Az[z])]<d-int(costs[z]) or lam[int(Az[z])]<minimum-int(costs[z]) for z in range(1,V)):
        return dict(status='EXCLUDED',reason='pure_Z')
    syn=np.array([sum(((form[a]&b).bit_count()%2)<<j for j,a in enumerate(A)) for b in range(N)])
    if len(set(map(int,syn)))!=V:raise ValueError('Nondegenerate dual pairing failed')
    W0=[int(np.flatnonzero(syn==(V-1 ^ (1<<i)))[0]) for i in range(p)]
    omega=[[((form[a]&b).bit_count()%2) for b in W0] for a in W0]
    # B=I+11^T, p even, B^-1=B. No independent diagonal gauges are removed.
    Bv=[v ^ ((V-1) if v.bit_count()%2 else 0) for v in range(V)]
    values=[np.array([w^int(Az[Bv[v]]) for v in range(V)]) for w in W0]
    for i in range(p):
        assert set(map(int,values[i]))==set(map(int,np.flatnonzero(syn==(V-1 ^ (1<<i)))))
    domains=[]
    for i in range(p):
        yes=allowed((1<<i)|(1<<(m-1)));domains.append(set(map(int,np.flatnonzero(yes[values[i]]))))
        if not domains[-1]:return dict(status='EXCLUDED',reason='unary',variable=i,W0=W0)
    relations={}
    for i in range(p):
        for j in range(i+1,p):
            yes=allowed((1<<i)|(1<<j));physical=yes[values[i][:,None]^values[j][None,:]]
            parity=np.array([[((vi>>j)^(vj>>i))&1 for vj in range(V)] for vi in range(V)])
            rel=physical & (parity==omega[i][j]);relations[i,j]=rel
            # Independent direct Gram comparison for every parametrized pair.
            direct=np.array([[(form[int(a)]&int(b)).bit_count()%2==0 for b in values[j]] for a in values[i]])
            assert np.array_equal(direct,parity==omega[i][j])
    changed=True;rounds=0
    while changed:
        changed=False;rounds+=1
        for (i,j),rel in relations.items():
            ai={v for v in domains[i] if any(rel[v,w] for w in domains[j])}
            aj={w for w in domains[j] if any(rel[v,w] for v in domains[i])}
            if not ai or not aj:return dict(status='EXCLUDED',reason='pair_arc_consistency',edge=[i,j],rounds=rounds)
            if ai!=domains[i] or aj!=domains[j]:domains[i],domains[j]=ai,aj;changed=True
    return dict(status='PASS_NEEDS_FULL_SEARCH',domain_sizes=list(map(len,domains)),W0=W0,
                domains=[sorted(v) for v in domains],rounds=rounds)


def main():
    inputs=json.loads((BASE.parent/'rank2_z_inputs.json').read_text())['classes'];reports=[];counts=Counter();start=time.monotonic()
    for item in inputs:
        idx=item['class_index'];data=json.loads(Path(item['table']).read_text())
        source=BASE/f'class_{idx}_z.jsonl'
        for lineno,line in enumerate(source.read_text().splitlines()):
            event=json.loads(line)
            if event['status']!='LEAF':continue
            result=screen(data,event['Z']);counts[result.get('reason',result['status'])]+=1
            record=dict(class_index=idx,line=lineno,Z=event['Z'],**result);reports.append(record)
            print(json.dumps(record),flush=True)
    result=dict(status='BOUNDED_LEAF_BATCH_FINISHED',tested=len(reports),counts=dict(counts),reports=reports,
        seconds=time.monotonic()-start,scope='Independent unary/pair necessary filters on saved bounded Z leaves; '
                    'passing leaves need a full solver, limited classes are not excluded.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE/'prepare.py']})
    (BASE/'pair_batch.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['reports','sha256']}))
if __name__=='__main__':main()
