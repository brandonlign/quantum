"""Complete native comparisons and genuine physical positive controls."""
import json,random,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,dual,gram,lift_coordinate,nullspace,reduced_basis,span,split,weight
from gauge_fix import supported_subspace
from fixed_suffix_search import quotient_input
from signed_extension import solve_extension
from rankone_light_budget.common import payload,OLD
from mixed_frames.stage3.common import atomic_json,sha
BINARY=HERE/'z_enum_projection'

def run(m,d,data,projection,binary=BINARY,xdom=False,limit=None):
    class Input:pass
    model=Input();model.m=m;model.d=d
    reps=[s for s,w in enumerate(data['distances']) if w>=d-1]
    cmd=[str(binary),'--gamma','--avoid-radical','--first-orbits']
    if projection:cmd+=['--projection']
    if xdom:cmd+=['--xdom']
    if limit is not None:cmd+=['--max-nodes',str(limit)]
    data=dict(data,labels=list(range(len(data['distances']))))
    p=subprocess.run(cmd,input=payload(model,data,reps),capture_output=True,text=True,check=True,timeout=30)
    st=json.loads(p.stderr);xs=[tuple(map(int,s.split()[1:])) for s in p.stdout.splitlines()]
    assert len(xs)==st['leaves']==len(set(xs));return st,xs

def valid(xs,G,m,d):
    full=list(xs)+[lift_coordinate((1<<len(xs))-1,xs)]
    return all(sum((c&x).bit_count()%2 for x in full)>=2*d-m and
               sum((c&x).bit_count()%2 for x in full)%2==0 for c in span(G) if c)

def main():
    records=[];rng=random.Random(17136)
    # Artificial complete trees test g=0..3 and different radical orientations.
    for b in range(3,7):
        for trial in range(4):
            pairs=(b-(trial%4))//2;g=b-2*pairs
            if g>3:continue
            basis=list(range(b));rng.shuffle(basis);Gamma=[0]*b
            for j in range(pairs):
                a,c=basis[2*j:2*j+2];Gamma[a]|=1<<c;Gamma[c]|=1<<a
            delta=[0]+[rng.randrange(2,5) for _ in range((1<<b)-1)]
            data=dict(distances=delta,Gamma=Gamma)
            G=nullspace(Gamma,b)
            for xdom in (False,True):
                old,L=run(4,3,data,False,binary=OLD,xdom=xdom)
                off,O=run(4,3,data,False,xdom=xdom)
                on,K=run(4,3,data,True,xdom=xdom)
                assert old['status']==off['status']==on['status']=='COMPLETE'
                assert L==O and old['nodes']==off['nodes']
                assert K==[x for x in L if valid(x,G,4,3)]
                records.append(dict(b=b,g=g,xdom=xdom,plain=len(L),retained=len(K),prunes=on['projection_pruned']))
    control8=json.loads((HERE/'positive_8.json').read_text())['result']['generators']
    # Higher distance thresholds force actual pruning, including g=3.
    for m,d,b,Gamma in [(4,4,3,[2,1,0]),(5,4,4,[2,1,0,0]),(5,4,5,[2,1,0,0,0])]:
        data=dict(distances=[0]+[d]*((1<<b)-1),Gamma=Gamma)
        G=nullspace(Gamma,b)
        old,L=run(m,d,data,False,binary=OLD)
        off,O=run(m,d,data,False)
        on,K=run(m,d,data,True)
        assert old['status']==off['status']==on['status']=='COMPLETE'
        assert L==O and old['nodes']==off['nodes'] and L
        assert K==[x for x in L if valid(x,G,m,d)] and len(K)<len(L)
        records.append(dict(m=m,d=d,b=b,g=len(G),plain=len(L),retained=len(K),prunes=on['projection_pruned']))
    models=[(ExtensionModel(5,0,3,3,[7<<3],[]),[]),(ExtensionModel(7,1,3,4,[15<<4],[7]),[])]
    C=[split(w,4,4)[1] for w in supported_subspace(dual(control8,8),8,range(4,8))]
    N=[split(w,4,4)[1] for w in supported_subspace(control8,8,range(4,8))]
    models.append((ExtensionModel(8,2,3,4,[15<<4],C),N))
    controls=[]
    for model,N in models:
        data=quotient_input(model,N);G=nullspace(data['Gamma'],len(data['Gamma']))
        st,L=run(model.m,model.d,data,False,binary=OLD)
        st,K=run(model.m,model.d,data,True)
        assert st['status']=='COMPLETE' and K==[x for x in L if valid(x,G,model.m,model.d)]
        positive=0
        for xs in L:
            answer=solve_extension(model,[lift_coordinate(x,data['units']) for x in xs],minimum_stabilizer_weight=model.d)
            if answer['status']!='FOUND':continue
            S,D=model.construct(answer['model_columns'])
            assert min(weight(w,model.n) for w in span(D) if w)>=model.d
            assert xs in K;positive+=1
        assert positive
        controls.append(dict(parameters=[model.n,model.k,model.d],g=len(G),positive_maps=positive,plain=len(L),retained=len(K)))
    assert [c['g'] for c in controls]==[0,1,2]
    for limit in (0,1):assert run(4,3,data,True,limit=limit)[0]['status']=='LIMIT'
    report=dict(status='PASS',complete_comparisons=records,controls=controls,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [BINARY,HERE/'z_enum_projection.c',Path(__file__),HERE/'PROOF.md',OLD,ROOT/'model.py',ROOT/'fixed_suffix_search.py',ROOT/'signed_extension.py']})
    atomic_json(HERE/'verification.json',report)
    print(json.dumps(dict(status='PASS',comparisons=len(records),controls=controls)),flush=True)
if __name__=='__main__':main()
