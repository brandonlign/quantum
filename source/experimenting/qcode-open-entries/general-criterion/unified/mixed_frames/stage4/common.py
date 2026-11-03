"""Stage-4-owned I/O and frozen, read-only stage-3 mathematical helpers."""
import json
from pathlib import Path
from mixed_frames.stage3.common import (ROOT, SOURCE, atomic_json, budget, case_input,
    graph_counts, input_text, physical_translate_costs, profile_mu, sha, target)
from .symmetry import automorphisms, coverage
from model import coordinates, lift_coordinate
from mixed_frames.stage2.common import graph_data
HERE=Path(__file__).resolve().parent
BINARY=HERE/'enumerate_z'
PROFILES={'888':[(0,1),(0,2),(1,0),(1,1),(2,0)]}


def hashes():
    paths=[SOURCE,ROOT/'model.py',*HERE.parent.glob('*.py'),
        * (HERE.parent/'stage2').glob('*.py'),* (HERE.parent/'stage3').glob('*.py'),
        HERE.parent/'stage3/enumerate_z',HERE.parent/'stage3/enumerate_z.cpp',
        HERE.parent/'stage3/manifest.json',HERE.parent/'stage3/overlap.json',
        HERE.parent/'stage3/PROOF.md',ROOT/'sol_workers/F_stage3_audit.md',
        *HERE.glob('*.py'),HERE/'enumerate_z.cpp',BINARY]
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}


def eligible(m,f,mu):
    labels=[s for s in range(1,1<<m.b) if m.distances[s]+f.costs[1]>=m.d]
    if not m.hull:
        weights={lift_coordinate(a,m.Gamma):w for a,w in enumerate(graph_data(m))}
        labels=[s for s in labels if weights[s]+f.costs[1]>=mu]
    return labels


def prepare():
    (HERE/'inputs').mkdir(exist_ok=True);(HERE/'certificates').mkdir(exist_ok=True)
    old=json.loads((HERE.parent/'stage3/manifest.json').read_text())
    overlap=json.loads((HERE.parent/'stage3/overlap.json').read_text())
    replay=json.loads((HERE/'provenance.json').read_text())
    assert replay['status']=='PASS' and replay['frozen_stage3_overlap_sha256']==sha(HERE.parent/'stage3/overlap.json')
    assert replay['replay_sha256']==sha(HERE/'overlap.json')
    assert json.loads((HERE/'overlap.json').read_text())==overlap
    assert overlap['status']=='PASS'
    assert sorted(r['class_index'] for r in overlap['mu7_reproduced'])==[7,32,41,50,52]
    assert len(overlap['mu8_reproduced'])==15 and len(overlap['singular_radicals'])==41
    provenance={p:sha(ROOT/p) for p in ['mixed_frames/stage3/overlap.json','mixed_frames/stage3/overlap.py','mixed_frames/stage3/PROOF.md','sol_workers/F_stage3_audit.md','mixed_frames/stage4/provenance.json','mixed_frames/stage4/overlap.json']}
    classes=[];skipped=[];autos=[]
    for index in range(56):
        m,f=target('888',index);aut=automorphisms(m)
        aut['class_index']=index;aut['P']=m.P;aut['C']=m.C;aut['Gamma']=m.Gamma
        atomic_json(HERE/f'certificates/aut_{index}.json',aut);autos.append(aut)
        maps=[g['syndrome_columns'] for g in aut['generators']]
        inverses=[[coordinates(1<<j,S) for j in range(m.b)] for S in maps]
        for profile in [(0,0),*PROFILES['888']]:
            key=f'888_{index}_{profile[0]}{profile[1]}'
            if profile==(0,0) or (profile[0]==0 and (m.hull or index in [7,32,41,50,52])):
                reason=('mu8 all hull-zero plus singular obstruction' if profile==(0,0) else
                    'mu7 singular balanced-kernel/radical/Plotkin obstruction' if m.hull else 'mu7 COMPLETE legacy zero-leaf replay')
                skipped.append(dict(key=key,class_index=index,profile=list(profile),reason=reason,provenance=provenance));continue
            row=dict(next(r for r in old['classes'] if r['key']==key))
            labels=eligible(m,f,row['mu_out']);orbits=coverage(labels,maps)
            cert=dict(key=key,class_index=index,profile=list(profile),eligible=labels,
                first_representatives=[o['representative'] for o in orbits],orbits=orbits,
                inverse_syndrome_columns=inverses,automorphism_file=f'certificates/aut_{index}.json',
                automorphism_sha256=sha(HERE/f'certificates/aut_{index}.json'),
                convention='T coefficient columns; S=T^{-T} syndrome columns; BFS forward witnesses; inverse columns for transport')
            file=HERE/f'certificates/{key}.json';atomic_json(file,cert)
            mask=HERE/f'inputs/{key}.sym';reps=set(cert['first_representatives'])
            mask.write_text(str(1<<m.b)+'\n'+' '.join(str(int(v in reps)) for v in range(1<<m.b))+'\n')
            path=HERE/row['input'];path.write_text(case_input(m,f,profile))
            assert sha(path)==row['input_sha256']
            row.update(first_images=cert['first_representatives'],symmetry=str(mask.relative_to(HERE)),
                symmetry_sha256=sha(mask),certificate=str(file.relative_to(HERE)),certificate_sha256=sha(file),
                eligible_labels=len(labels),canonical=True)
            classes.append(row)
        print(json.dumps(dict(class_index=index,hull=m.hull,generators=len(maps),order=aut['graph_group_order'])),flush=True)
    atomic_json(HERE/'exclusions.json',dict(status='INHERITED',cases=skipped,provenance=provenance,scope='No new DFS exclusion; stage-3 evidence trusted at recorded hashes.'))
    atomic_json(HERE/'manifest.json',dict(classes=classes,hashes=hashes(),
        exclusions_sha256=sha(HERE/'exclusions.json'),
        partitions=sum(len(r['first_images']) for r in classes),
        scope='888 only; combined lex-minimum necessary conditions; serial explicit launch'))
    return dict(cases=len(classes),skipped=len(skipped),partitions=sum(len(r['first_images']) for r in classes))


if __name__=='__main__':print(json.dumps(prepare()))
