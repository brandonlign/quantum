"""Freeze distinguished-prefix and suffix root coverage for the state bases."""
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from suffix_symmetries import automorphisms
from gauge_rankone_h5.common import (PREFIX,SOURCE,prefix_candidates,cheap_frame,
    target,quotient,input_text,suffix_generators,first_orbits,atomic_json,sha)
from model import span,weight


def main():
    parent=json.loads(PREFIX.read_text())['groups']['4']['classes']
    all_prefix=prefix_candidates(parent)
    selected=[];coverage=[]
    for index,row in enumerate(parent):
        H=row['generators'];maps,order=automorphisms(H,8)
        remaining={w for w in span(H) if weight(w,8)==8};orbits=[]
        while remaining:
            root=min(remaining);seen={root};todo=[root]
            while todo:
                v=todo.pop()
                for g in maps:
                    w=g(v);assert w in remaining
                    if w not in seen:seen.add(w);todo.append(w)
            remaining-=seen
            pr=next(p for p in all_prefix if p['class_index']==index and p['distinguished_word']==root)
            selected.append(pr);orbits.append(dict(prefix_index=pr['prefix_index'],representative=root,members=sorted(seen)))
        coverage.append(dict(class_index=index,physical_group_order=order,orbits=orbits))
    assert len(selected)==7
    (HERE/'inputs').mkdir(exist_ok=True);cases=[];suffixes=[]
    source=json.loads(SOURCE.read_text())['classes']
    frames={p['prefix_index']:cheap_frame(8,p['R'],f'h5_gauge_{p["prefix_index"]}') for p in selected}
    for c in source:
        if c['representative']['hull']!=5:continue
        C=c['representative']['generators'];p0=selected[0]
        model,_=target(tuple(C),tuple(p0['R']),f'h5_gauge_{p0["prefix_index"]}')
        data=quotient(model);generators,order=suffix_generators(model,data)
        suffixes.append(dict(cls=c['class_index'],C=model.C,data=data,syndrome_generators=generators,physical_group_order=order))
        for pr in selected:
            f=frames[pr['prefix_index']]
            key=f'h5_{c["class_index"]}_p{pr["prefix_index"]}'
            path=HERE/'inputs'/f'{key}.txt';path.write_text(input_text(model,f,data))
            eligible=[s for s in range(256) if s and data['distances'][s]+f.costs[1]>=6]
            orbits=first_orbits(eligible,generators)
            cases.append(dict(key=key,cls=c['class_index'],prefix_index=pr['prefix_index'],C=model.C,R=f.R,E=f.E,
                              input=str(path.relative_to(HERE)),input_sha256=sha(path),
                              first_roots=[o['representative'] for o in orbits],eligible=eligible,orbits=orbits,
                              prefix_group=f'prefix_masks/{pr["prefix_index"]}.txt',
                              prefix_group_sha256=sha(HERE/'prefix_masks'/f'{pr["prefix_index"]}.txt')))
    assert len(suffixes)==118 and len(cases)==826
    atomic_json(HERE/'suffix_quotients.json',suffixes)
    atomic_json(HERE/'prefix_word_coverage.json',coverage)
    result=dict(status='PREPARED',cases=cases,prefixes=selected,classes=118,
                partitions=sum(len(c['first_roots']) for c in cases),
                source_sha256=sha(SOURCE),prefix_source_sha256=sha(PREFIX),
                scope='Pure state bases of the rank-one hull-five column; inverse selection is mandatory')
    atomic_json(HERE/'manifest.json',result)
    print(json.dumps(dict(status='PREPARED',prefix_word_orbits=7,cases=826,partitions=result['partitions'])),flush=True)


if __name__=='__main__':main()
