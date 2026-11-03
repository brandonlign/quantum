"""Read-only legacy evidence replay and independent physical table agreement.

This reproduces inherited exclusions, not a new exhaustive stage-3 rerun.
"""
from collections import Counter
import json
from model import coordinates,gram,lift_coordinate,weight,reduced_basis,nullspace
from mixed_frames.frames import xor_table
from mixed_frames.suffix_tables import verify_physical
from .common import HERE,ROOT,atomic_json,sha,target


def check_hashes(d):
    for path,digest in d['sha256'].items():assert sha(path)==digest,(path,'changed legacy evidence')

def orbit_reps(high,maps):
    remaining=set(high);required=[]
    while remaining:
        seed=min(remaining);orbit={seed};todo=[seed]
        while todo:
            v=todo.pop()
            for T in maps:
                w=lift_coordinate(v,T);assert w in high
                if w not in orbit:orbit.add(w);todo.append(w)
        assert orbit<=remaining;remaining-=orbit;required.append(seed)
    return required

def flags(delta,lam,mu):
    return [int(a!=0)*(int(d>=4 and w>=mu-2)+2*int(d>=2 and w>=mu-4)) for a,(d,w) in enumerate(zip(delta,lam))]

def mask(delta,lam,mu):
    return [sum(1<<c for c in range(9) if d+c>=6 and w+c>=mu) for d,w in zip(delta,lam)]

def check_input(path,table,mu,required=None):
    lines=path.read_text().splitlines();header=list(map(int,lines[0].split()))
    assert header[:2]==[4096,6] and list(map(int,lines[1].split()))==[1,1,2,2,1,1]
    fs=flags(table['coset_distances'],table['suffix_weights'],mu)
    masks=mask(table['coset_distances'],table['suffix_weights'],mu)
    for a,line in enumerate(lines[2:4098]):
        row=list(map(int,line.split()))
        assert row[:2]==[fs[a],table['form_rows'][a]]
        if len(row)>2:assert row[2]==masks[a]
    if required is not None:
        assert int(lines[4098])==len(required) and list(map(int,lines[4099].split()))==required

def main():
    auto=json.loads((ROOT/'rank2_z/aut_orbits.json').read_text())
    tables={};mapdata={r['class_index']:r for r in auto['classes']}
    result=dict(status='PASS',hull_zero=[],physical_suffix_checks=[],singular_radicals=[],mu7_reproduced=[],mu8_reproduced=[])
    saved=json.loads((HERE/'suffix_tables.json').read_text())['classes']
    for t in saved:
        # Independent full physical fiber enumeration, all 56 hull regimes.
        check=verify_physical(t)
        result['physical_suffix_checks'].append(dict(class_index=t['class_index'],**check))
        m,f=target('888',t['class_index'])
        assert t['P']==m.P and t['Gamma']==m.Gamma
        if m.hull:
            radical=[lift_coordinate(a,m.P) for a in xor_table(nullspace(m.Gamma,m.b))]
            weights=Counter(weight(w,9) for w in radical if w)
            assert len(radical)==1<<m.hull
            if m.hull==2:assert weights=={6:3}
            else:assert m.hull==6 and set(weights)<={6,8}
            result['singular_radicals'].append(dict(class_index=t['class_index'],hull=m.hull,nonzero_weights=dict(weights)))
            continue
        index=t['class_index'];table=json.loads((ROOT/f'rank2_z/class_{index}.json').read_text())
        words=xor_table(m.P);lam=[weight(w,9) for w in words]
        forms=[m.syndrome(w) for w in words]
        delta=[m.distances[s] for s in forms]
        assert table['P_generators']==m.P and table['C_generators']==m.C
        assert table['form_rows']==forms and table['suffix_weights']==lam and table['coset_distances']==delta
        assert table['flags']==flags(delta,lam,8)
        assert table['needs']==[1,1,2,2,1,1]
        # Independent physical old basis and all 4096 prefix class costs.
        oldZ=[((1<<i)|(1<<7))<<8 for i in range(6)]
        oldF=[(1<<i)|(1<<7) for i in range(6)]
        basis=oldZ+oldF
        assert gram(basis,8)==[(1<<(i+6))^(63<<6) for i in range(6)]+[(1<<i)^63 for i in range(6)]
        change=[coordinates(v,f.R+f.E)>>2 for v in basis]
        for e,w in enumerate(xor_table(basis)):
            wt=min(weight(w^r,8) for r in xor_table(f.R))
            assert f.costs[lift_coordinate(e,change)]==wt
        # Recheck each linear table symmetry; it preserves the full criterion.
        sym=mapdata[index];assert sym['table_sha256']==sha(ROOT/f'rank2_z/class_{index}.json')
        maps=sym['verified_generator_coefficient_maps']
        for T in maps:
            assert len(reduced_basis(T))==12
            assert gram([words[v] for v in T],9)==m.Gamma
            for a,ta in enumerate(xor_table(T)):
                assert lam[a]==lam[ta] and delta[a]==delta[ta]
        tables[index]=table
        result['hull_zero'].append(dict(class_index=index,mu8_high=sum(bool(x&1) for x in table['flags']),
             mu7_high=sum(bool(x&1) for x in flags(delta,lam,7)),prefix_classes_checked=4096,verified_maps=len(maps)))
    for index in [7,32,41,50,52]:
        table=tables[index];high=[a for a,(d,w) in enumerate(zip(table['coset_distances'],table['suffix_weights'])) if a and d>=4 and w>=5]
        required=orbit_reps(high,mapdata[index]['verified_generator_coefficient_maps'])
        if index in (32,52):
            source=ROOT/f'rank2_mu7_orbits/class_{index}.txt'
            screen=json.loads((ROOT/'rank2_mu7_orbits/screen.json').read_text())
            r=next(r for r in screen['records'] if r['class_index']==index)
            assert r['status']=='COMPLETE' and r['Z_excluded'] and r['leaves']==0
            assert sha(source)==r['input_sha256'] and sha(ROOT/f'rank2_z/class_{index}.json')==r['source_sha256']
            check_hashes(screen)
            check_input(source,table,7,required)
            events=[json.loads(l) for l in Path(r['log']).read_text().splitlines()]
            assert len(events)==1 and events[-1]['status']=='COMPLETE' and events[-1]['leaves']==0
            assert events[-1]['nodes']==r['nodes']
            result['mu7_reproduced'].append(dict(class_index=index,profiles=[[0,1],[0,2]],type='legacy_complete_zero_leaf',first_orbits=len(required),nodes=r['nodes'],input_sha256=sha(source),log_sha256=sha(r['log'])))
        else:
            suffix='physical_partitions' if index in (7,41) else 'partitions'
            path=ROOT/f'rank2_mu7/class{index}_{suffix}.json'
            d=json.loads(path.read_text());check_hashes(d)
            assert d['status']=='EXCLUDED_MIN7_CANDIDATE' and d['minimum_stabilizer']==7 and d['required_first']==required
            covered=set();nodes=0
            for r in d['records']:
                assert r['status']=='PARTITION_EXCLUDED' and r['found'] is None
                assert r['z']['status']=='COMPLETE' and r['z']['leaves']==0 and not r['verdicts']
                assert not covered.intersection(r['first']) and len(r['first'])==len(set(r['first']))
                covered.update(r['first']);nodes+=r['z']['nodes']
                assert sha(r['input'])==r['input_sha256'] and sha(r['log'])==r['log_sha256']
                check_input(Path(r['input']),table,7,r['first'])
                events=[json.loads(l) for l in Path(r['log']).read_text().splitlines()]
                assert events==[r['z']]
            assert sorted(covered)==required==d['covered_first']
            result['mu7_reproduced'].append(dict(class_index=index,profiles=[[0,1],[0,2]],type='legacy_complete_partition_replay',partitions=len(d['records']),first_orbits=len(required),nodes=nodes,evidence_sha256=sha(path)))
    branch=json.loads((ROOT/'rank2_branch_result.json').read_text());check_hashes(branch)
    assert branch['status']=='COMPLETE_BRANCH_EXCLUSION' and len(branch['records'])==15
    for name in ('symmetric_long','remaining_long'):
        batch=json.loads((ROOT/f'rank2_z/{name}.json').read_text())
        assert batch['symmetry_first_orbits']
        for file,digest in batch['sha256'].items():assert sha(ROOT/'rank2_z'/file)==digest
    for r in branch['records']:
        index=r['class_index'];table=tables[index]
        assert r['status']=='COMPLETE' and r['Z_excluded'] and r['leaves']==0
        log=Path(r['log']);source=log.with_suffix('.txt')
        assert sha(log)==r['log_sha256'] and sha(source)==r['input_sha256']
        required=orbit_reps([a for a,flag in enumerate(table['flags']) if flag&1],mapdata[index]['verified_generator_coefficient_maps'])
        check_input(source,table,8,required)
        events=[json.loads(l) for l in log.read_text().splitlines()]
        assert len(events)==1 and events[-1]['status']=='COMPLETE' and events[-1]['leaves']==0 and events[-1]['nodes']==r['nodes']
        result['mu8_reproduced'].append(dict(class_index=index,profile=[0,0],L=0,nodes=r['nodes'],input_sha256=sha(source),log_sha256=sha(log)))
    assert len(result['hull_zero'])==15
    result['scope']='Independent physical tables/basis/linear table symmetries and saved complete evidence replay. Inherited mu7/8 exclusions reproduce exact-profile consequences; no new exhaustive DFS replay.'
    result['singular_mu7_mu8']='hull2 radicals are all wt6, so no >=7 one-dimensional N; hull6 requires dim N=3 and |N|=8, exceeding the quaternary d8 length9 Plotkin bound floor(32/5)=6. Only hull0 remains for profiles with A6=0.'
    atomic_json(HERE/'overlap.json',result)
    print(json.dumps(dict(status='PASS',physical_classes=56,hull_zero=15,mu7=[r['class_index'] for r in result['mu7_reproduced']],mu8=len(result['mu8_reproduced']))))

from pathlib import Path
if __name__=='__main__':main()
