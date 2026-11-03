"""Prepare exact rank-two Z-stage tables for target17's15 hull-zero classes."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import ExtensionModel,span,lift_coordinate,weight
BASE=Path(__file__).resolve().parent.parent
OUT=BASE/'rank2_z'


def table(C,m,n,k,d,minimum):
    R=[(1<<m)-1,((1<<m)-1)<<m]
    model=ExtensionModel(n,k,d,m,R,C)
    p=m-2
    if model.hull!=0 or model.b!=2*p:raise ValueError('Requires nondegenerate equal-size prefix and suffix quotient')
    # Coefficient labels in P, not syndrome labels. Syndrome(p_alpha)=Gamma alpha.
    words=[lift_coordinate(a,model.P) for a in range(1<<model.b)]
    forms=[model.syndrome(w) for w in words]
    flags=[]
    for alpha,w in enumerate(words):
        delta=model.distances[forms[alpha]];lam=weight(w,model.f)
        flag=int(delta>=d-2 and lam>=minimum-2) | (2*int(delta>=d-4 and lam>=minimum-4))
        flags.append(flag if alpha else 0)
    needs=[]
    for j in range(1,p+1):
        support=j+(j%2);cost=min(support,m-support)
        if cost not in [2,4]:raise ValueError('Backend supports pure-Z costs2/4')
        needs.append(1 if cost==2 else 2)
    return model,dict(flags=flags,form_rows=forms,needs=needs,
                      P_generators=model.P,C_generators=model.C,
                      suffix_weights=[weight(w,model.f) for w in words],
                      coset_distances=[model.distances[s] for s in forms])


def write_input(path,data,p,cap,millis,leafcap):
    rows=[f'{len(data["flags"])} {p} {cap} {millis} {leafcap}', ' '.join(map(str,data['needs']))]
    rows += [f'{f} {g}' for f,g in zip(data['flags'],data['form_rows'])]
    path.write_text('\n'.join(rows)+'\n')


def main():
    source=BASE/'suffix9_classes.json';census=json.loads(source.read_text());items=[]
    for v in census['classes']:
        if v['representative']['hull']!=0:continue
        idx=v['class_index'];model,data=table(v['representative']['generators'],8,17,3,6,8)
        data.update(class_index=idx,parameters=[17,3,6],m=8,r=2,p=6,
                    input_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
        path=OUT/f'class_{idx}.json';path.write_text(json.dumps(data)+'\n')
        write_input(OUT/f'class_{idx}.txt',data,6,1000000,3000,20)
        items.append(dict(class_index=idx,table=str(path),input=str(OUT/f'class_{idx}.txt'),
                          high_labels=sum(bool(f&1) for f in data['flags'])))
    assert len(items)==15
    (BASE/'rank2_z_inputs.json').write_text(json.dumps(dict(status='PREPARED',classes=items,
        scope='Necessary rank-two pure-Z constraints, including minimum stabilizer weight8; not full extensions.'),indent=2)+'\n')
    print(json.dumps(items))
if __name__=='__main__':main()
