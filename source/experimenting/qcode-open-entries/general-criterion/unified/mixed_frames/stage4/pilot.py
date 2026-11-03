"""Handful of COMPLETE small frontiers; limits prevent any exhaustive batch."""
import json
import time
from mixed_frames.stage3.common import BINARY as OLD
from mixed_frames.stage3.frontier_pilot import legacy_frontiers
from .common import HERE,atomic_json,target
from .checks import cpp
from .symmetry import joint_normal


def main():
    records=[];generations=[];started=time.monotonic()
    seen=set();frontiers=[]
    for item in legacy_frontiers():
        i=item['class_index'];key=f'888_{i}_20';cert=json.loads((HERE/f'certificates/{key}.json').read_text())
        normal=joint_normal(item['X'],cert)
        if (i,tuple(normal)) in seen:continue
        seen.add((i,tuple(normal)));frontiers.append((i,normal[:5],(2,0),'legacy leaf canonical depth-five'))
        if len(frontiers)==6:break
    # A few larger complete subtrees, plus one singular case in each hull regime.
    for i in [25,34]:
        normal=next(p for j,p,_,_ in frontiers if i==j)
        for depth in [4,3,2]:
            frontiers.append((i,normal[:depth],(2,0),f'canonical depth-{depth}'))
    for i in [0,12]:
        m,f=target('888',i);key=f'888_{i}_20'
        r=cpp(m,f,mu=6,counts=(2,0),depth=5,nodes=3000,seconds=.5,leaves=1,symmetry=HERE/f'inputs/{key}.sym')
        generations.append(dict(class_index=i,**r))
        if r['leaves']:
            for depth in [5,4,3]:frontiers.append((i,r['leaves'][0][:depth],(2,0),f'singular generated depth-{depth}'))
    for i,prefix,counts,origin in frontiers:
        m,f=target('888',i);key=f'888_{i}_{counts[0]}{counts[1]}'
        pair={}
        for lane,binary in [('stage3',OLD),('stage4',None)]:
            kw=dict(mu=6,counts=counts,prefix=prefix,nodes=30000,seconds=1.5)
            if binary is not None:kw['binary']=binary
            else:kw['symmetry']=HERE/f'inputs/{key}.sym'
            r=cpp(m,f,**kw);pair[lane]=r
        complete=all(r['complete'] for r in pair.values())
        row=dict(class_index=i,hull=m.hull,profile=list(counts),prefix=prefix,origin=origin,
            both_complete=complete,**pair)
        records.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ['stage3','stage4']}|{k:{q:v for q,v in r.items() if q!='leaves'} for k,r in pair.items()}),flush=True)
    complete=[r for r in records if r['both_complete']]
    assert len(complete)>=5
    totals={lane:{k:sum(r[lane][k] for r in complete) for k in ['nodes','candidates','leaf_count','seconds']} for lane in ['stage3','stage4']}
    atomic_json(HERE/'pilots.json',dict(records=records,generation=generations,complete_frontiers=len(complete),
        totals=totals,seconds=time.monotonic()-started,
        scope='Identical specified canonical prefixes; COMPLETE remaining subtrees only. Not whole root partitions or classes; no Y completion. Nonrandom survival-selected sample.'))
    print(json.dumps(dict(complete_frontiers=len(complete),totals=totals)),flush=True)


if __name__=='__main__':main()
