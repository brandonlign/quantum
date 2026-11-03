"""Frozen stage-3 physical tables and exact profile predicates."""
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from model import ExtensionModel, lift_coordinate, weight, span
from mixed_frames.frames import xor_table
from mixed_frames.suffix_tables import SOURCE, build
from mixed_frames.stage2.common import graph_data, graph_minimum, graph_z_ok, graph_partial_ok
from .frames import mixed_frame

if not __debug__:
    raise RuntimeError('Stage 3 requires assertions; do not use Python -O')

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BINARY=HERE/'enumerate_z'
PROFILES={'688':[(1,0),(1,1),(2,0)], '888':[(0,1),(0,2),(1,0),(1,1),(2,0)]}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def hashes():
    paths=[SOURCE,ROOT/'model.py',*HERE.parent.glob('*.py'),
           HERE.parent/'suffix_tables.json',ROOT/'weight8_split/frame_count_certificates.json',
           ROOT/'mixed_frames/stage2/common.py', *HERE.glob('*.py'),HERE/'enumerate_z.cpp',BINARY]
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}

def atomic_json(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)

@lru_cache(None)
def target(name,index):
    frame=mixed_frame(name)
    c=next(c for c in json.loads(SOURCE.read_text())['classes'] if c['class_index']==index)
    return ExtensionModel(17,3,6,8,frame.R,c['representative']['generators']),frame

@lru_cache(None)
def translate_costs(name):
    frame=mixed_frame(name)
    return physical_translate_costs(frame)

def physical_translate_costs(frame):
    rs=xor_table(frame.R)
    return [[weight(w^r,frame.m) for r in rs] for w in xor_table(frame.E)]

def budget(frame,profile):
    r=Counter(weight(v,frame.m) for v in xor_table(frame.R)[1:])
    b6,b7=profile[0]-r[6],profile[1]-r[7]
    assert b6>=0 and b7>=0
    return b6,b7

def graph_counts(model,frame,columns,w6=6,w7=7):
    pre=physical_translate_costs(frame);ws=graph_data(model);p=frame.p
    hist=Counter()
    for a in range(1,1<<model.b):
        e=sum(((a&c).bit_count()&1)<<((i+p)%(2*p)) for i,c in enumerate(columns))
        hist.update(t+ws[a] for t in pre[e])
    return dict(minimum=min(hist,default=model.n+1),light=sum(v for w,v in hist.items() if w<=w7),
                n6=hist[w6],n7=hist[w7],histogram={str(w):n for w,n in sorted(hist.items())})

def input_text(model,frame,mu=0,rank=True,domains=True,counts=None,w6=6,w7=7):
    arrays=[[frame.p,model.b,model.d,mu,int(rank),int(domains)],frame.costs,
            model.distances,model.Gamma,graph_data(model)]
    if counts is not None:
        b6,b7=counts
        assert mu>=w6
        pre=physical_translate_costs(frame)
        arrays += [[b6+b7,b6,b7,w7+1,w6,w7,len(pre[0])], [v for row in pre for v in row]]
    return '\n'.join(' '.join(map(str,a)) for a in arrays)+'\n'

def profile_mu(frame,profile):
    b6,b7=budget(frame,profile)
    return 6 if b6 else (7 if b7 else 8)

def case_input(model,frame,profile):
    return input_text(model,frame,profile_mu(frame,profile),counts=budget(frame,profile))

def prefilter(model,frame):
    # Injective X maps distinct nonzero pure-Z classes to distinct syndromes.
    assert max(frame.cost[0])<model.d
    rows=[]
    for c in sorted(set(frame.cost[0][1:])):
        required=sum(0<v<=c for v in frame.cost[0][1:])
        available=sum(s!=0 and delta>=model.d-c for s,delta in enumerate(model.distances))
        rows.append(dict(cost_at_most=c,required=required,available=available))
    return dict(tests=rows,verdict='PREFILTER_REJECT' if any(r['available']<r['required'] for r in rows) else 'CANDIDATE')

def prepare():
    (HERE/'inputs').mkdir(exist_ok=True)
    tables=[];records=[];classes=[]
    saved=json.loads((HERE.parent/'suffix_tables.json').read_text())['classes']
    frame_records=[]
    for name in PROFILES:
        frame=mixed_frame(name);frame.verify_physical()
        frame_records.append(dict(**frame.record(),translate_weights=physical_translate_costs(frame)))
    for index in range(56):
        model,frame=target('888',index)
        fresh=build(model.C)
        assert all(saved[index][k]==v for k,v in fresh.items())
        tables.append(dict(class_index=index,**fresh,coefficient_weights=graph_data(model)))
        for name,profiles in PROFILES.items():
            m,f=target(name,index);cap=prefilter(m,f)
            records.append(dict(frame=name,class_index=index,**cap))
            if cap['verdict']!='CANDIDATE':continue
            for profile in profiles:
                key=f'{name}_{index}_{profile[0]}{profile[1]}'
                path=HERE/'inputs'/(key+'.txt');path.write_text(case_input(m,f,profile))
                classes.append(dict(key=key,frame=name,class_index=index,profile=list(profile),
                    budget=list(budget(f,profile)),d=6,mu_out=profile_mu(f,profile),
                    input=str(path.relative_to(HERE)),input_sha256=sha(path),
                    first_images=[s for s in range(1<<m.b) if s and m.distances[s]+f.costs[1]>=6],
                    R=f.R,E=f.E,C=m.C,P=m.P,Gamma=m.Gamma,hull=m.hull))
    atomic_json(HERE/'prefix_tables.json',dict(frames=frame_records))
    atomic_json(HERE/'suffix_tables.json',dict(source_sha256=sha(SOURCE),classes=tables))
    atomic_json(HERE/'capacity_certificate.json',dict(records=records,
        proof='Nonzero pure-Z costs <6 force injectivity. For every c, #classes of cost <=c cannot exceed #nonzero syndromes of distance >=6-c.'))
    atomic_json(HERE/'manifest.json',dict(classes=classes,hashes=hashes(),
        profiles={k:[list(p) for p in v] for k,v in PROFILES.items()},
        table_sha256={p.name:sha(p) for p in [HERE/'prefix_tables.json',HERE/'suffix_tables.json',HERE/'capacity_certificate.json']}))
    return dict(cases=len(classes),capacity_rejected={name:[r['class_index'] for r in records if r['frame']==name and r['verdict']=='PREFILTER_REJECT'] for name in PROFILES})

if __name__=='__main__':print(json.dumps(prepare()))
