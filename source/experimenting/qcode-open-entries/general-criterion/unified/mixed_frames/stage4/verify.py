"""Bounded, serial orbit coverage tests against the frozen Stage-3 binary."""
import itertools
import json
import time
import tempfile
import copy
from . import run
from model import ExtensionModel,coordinates,lift_coordinate,span,symp
from mixed_frames.frames import adapted_frame,xor_table
from mixed_frames.stage3.common import BINARY as OLD
from mixed_frames.stage3.frontier_pilot import legacy_frontiers
from .common import HERE,BINARY,atomic_json,budget,eligible,graph_counts,hashes,sha,target
from .checks import cpp
from .symmetry import (affine_points,anchor_normal,automorphisms,coverage,joint_normal,
    physical_apply,validate_certificate)
from .run import complete_leaf


def fixture_cert(m,f,mu,name):
    aut=automorphisms(m);maps=[g['syndrome_columns'] for g in aut['generators']]
    labs=eligible(m,f,mu);os=coverage(labs,maps)
    c=dict(eligible=labs,orbits=os,first_representatives=[o['representative'] for o in os],
        inverse_syndrome_columns=[[coordinates(1<<i,S) for i in range(m.b)] for S in maps])
    validate_certificate(c,aut)
    mask=HERE/f'{name}.sym';reps=set(c['first_representatives'])
    mask.write_text(str(1<<m.b)+'\n'+' '.join(str(int(a in reps)) for a in range(1<<m.b))+'\n')
    atomic_json(HERE/f'{name}.certificate.json',dict(certificate=c,automorphisms=aut))
    return c,mask


def prefix_check():
    m,f=target('888',7);weights=__import__('mixed_frames.stage3.common',fromlist=['physical_translate_costs']).physical_translate_costs(f)
    star=[coordinates(((1<<i)|(1<<7))<<8,f.R+f.E)>>2 for i in range(7)]
    assert star==[1,2,4,8,16,32,63]
    result=[]
    for k in range(8):
        sites=list(range(8));letters=[[0,1,2,3] for _ in sites]
        if k<7:sites[k],sites[k+1]=sites[k+1],sites[k]
        else:letters=[[0,3,2,1] for _ in sites] # X<->Y, Z fixed
        apply=lambda w:physical_apply(w,sites,letters,8)
        assert {apply(w) for w in xor_table(f.R)}==set(xor_table(f.R))
        U=[coordinates(apply(w),f.R+f.E)>>2 for w in f.E]
        assert all(symp(apply(a),apply(b),8)==symp(a,b,8) for a in f.E for b in f.E)
        for e in range(4096):
            ue=lift_coordinate(e,U)
            assert f.costs[e]==f.costs[ue] and sorted(weights[e])==sorted(weights[ue])
        assert all(u<64 for u in U[:6])
        if k==7:assert U[:6]==[1<<i for i in range(6)]
        result.append(dict(generator=k,quotient_columns=U,classes_checked=4096))
    return result


def main():
    started=time.monotonic();report=dict(status='PASS',prefix_generators=prefix_check(),fixtures=[],target_corpus=[],bounded_prefixes=[],positive_controls=[])
    # Baseline disabled-symmetry lane must exactly preserve the predecessor.
    for m in [ExtensionModel(6,0,3,4,[15,240],[]),ExtensionModel(7,1,2,4,[240,15],[24,3])]:
        f=adapted_frame(4,m.R)
        for mu in [0,m.d,m.d+1]:
            cert,mask=fixture_cert(m,f,mu,f'fixture_{m.n}_{mu}')
            for rank,domains in itertools.product([False,True],repeat=2):
                kw=dict(mu=mu,rank=rank,domains=domains)
                a=cpp(m,f,binary=OLD,**kw);unreduced=cpp(m,f,canonical=False,**kw)
                assert (a['leaves'],a['nodes'],a['status'])==(unreduced['leaves'],unreduced['nodes'],unreduced['status'])
                b=cpp(m,f,symmetry=mask,**kw)
                assert a['complete'] and b['complete']
                leaves={tuple(xs) for xs in b['leaves']}
                normalized={tuple(joint_normal(xs,cert)) for xs in a['leaves']}
                assert normalized<=leaves
                report['fixtures'].append(dict(n=m.n,hull=m.hull,mu=mu,rank=rank,domains=domains,
                    stage3_nodes=a['nodes'],stage4_nodes=b['nodes'],stage3_leaves=a['leaf_count'],stage4_leaves=b['leaf_count'],covered=len(normalized)))
        # A genuine completed physical positive control, including exact counts.
        original=[3,12] if m.n==6 else [4,8]
        positives=complete_leaf(m,f,original,m.d,None,m.d,m.d+1)
        assert positives['accepted']
        cert,mask=fixture_cert(m,f,m.d,f'positive_{m.n}');normal=joint_normal(original,cert)
        completed=complete_leaf(m,f,normal,m.d,None,m.d,m.d+1)
        assert completed['accepted']
        for measure in positives['measurements']:
            if measure['minimum']<m.d:continue
            counts=(measure['n6'],measure['n7'])
            enumerated=cpp(m,f,mu=m.d,counts=counts,w6=m.d,w7=m.d+1,prefix=normal,symmetry=mask)
            assert enumerated['leaves']==[normal]
            matching=complete_leaf(m,f,normal,m.d,counts,m.d,m.d+1)
            assert matching['accepted']
        report['positive_controls'].append(dict(n=m.n,hull=m.hull,original_X=original,representative=normal,
            original_solutions=len(positives['accepted']),representative_solutions=len(completed['accepted']),exact_profiles_survive=True))
    # Exercise the actual Stage-4 streaming worker, canonical mask and OK ack.
    m=ExtensionModel(6,0,3,4,[15,240],[]);f=adapted_frame(4,m.R)
    cert,mask=fixture_cert(m,f,3,'worker_fixture')
    old_here,old_target=run.HERE,run.target
    try:
        with tempfile.TemporaryDirectory(dir=HERE,prefix='worker_') as directory:
            folder=__import__('pathlib').Path(directory);run.HERE=folder;run.target=lambda name,index:(m,f)
            (folder/'input.txt').write_text(__import__('mixed_frames.stage4.common',fromlist=['input_text']).input_text(m,f,3))
            (folder/'mask.sym').write_bytes(mask.read_bytes())
            row=dict(key='fixture',frame='fixture',class_index=0,profile=[0,0],budget=None,d=3,mu_out=3,
                input='input.txt',input_sha256=sha(folder/'input.txt'),symmetry='mask.sym',
                symmetry_sha256=sha(folder/'mask.sym'),certificate_sha256='fixture')
            result=run.partition((row,1,{'fixture':'stage4'}))
            assert result['status']=='COMPLETE' and result['surviving_solutions']>0 and not result['excluded'],result
            path=folder/'runs/fixture/1/checkpoint.json';attack=copy.deepcopy(result)
            attack['surviving_solutions']=0;attack['excluded']=True;atomic_json(path,attack)
            try:run.verify_partition(path.parent,row,1,{'fixture':'stage4'})
            except AssertionError:pass
            else:raise AssertionError('False zero-survivor checkpoint accepted')
            report['worker_control']=dict(status=result['status'],leaves=result['leaves'],
                survivors=result['surviving_solutions'],false_exclusion_rejected=True)
    finally:run.HERE,run.target=old_here,old_target
    # All 19 previously found actual 888 mu6 leaves, with unchanged and finite budgets.
    for item in legacy_frontiers():
        m,f=target('888',item['class_index']);key=f'888_{item["class_index"]}_20'
        cert=json.loads((HERE/f'certificates/{key}.json').read_text());normal=joint_normal(item['X'],cert)
        for profile in [None,(1,0),(1,1),(2,0)]:
            kw=dict(mu=6,counts=profile)
            a=cpp(m,f,prefix=item['X'],binary=OLD,**kw)
            b=cpp(m,f,prefix=normal,symmetry=HERE/f'inputs/{key}.sym',**kw)
            assert a['complete'] and b['complete']
            assert not a['leaves'] or b['leaves']==[normal]
            # Here full-profile sector counts are invariant, giving both directions.
            assert bool(a['leaves'])==bool(b['leaves'])
            report['target_corpus'].append(dict(class_index=item['class_index'],hull=m.hull,profile=profile,
                source_X=item['X'],representative=normal,source_sha256=item['source_sha256'],
                stage3_leaves=a['leaf_count'],stage4_leaves=b['leaf_count'],covered=True))
    # Search small depth-five target frontiers, then replay every observed full leaf
    # at its constructive orbit representative. LIMIT means sampled coverage only.
    for item in legacy_frontiers()[:5]:
        index=item['class_index'];m,f=target('888',index);key=f'888_{index}_20'
        cert=json.loads((HERE/f'certificates/{key}.json').read_text())
        a=cpp(m,f,mu=6,prefix=item['X'][:5],binary=OLD,nodes=2000,seconds=.4,leaves=12)
        reps=[]
        for xs in a['leaves']:
            normal=joint_normal(xs,cert)
            b=cpp(m,f,mu=6,prefix=normal,symmetry=HERE/f'inputs/{key}.sym',nodes=8,seconds=.3)
            assert b['complete'] and b['leaves']==[normal];reps.append(normal)
        report['bounded_prefixes'].append(dict(class_index=index,hull=m.hull,prefix=item['X'][:5],
            stage3={k:v for k,v in a.items() if k!='leaves'},observed_leaves=a['leaves'],representatives=reps,
            covered=len(reps),scope='Each observed leaf replayed with COMPLETE fixed-X test; LIMIT subtree remains unresolved'))
    # Short root samples explicitly include all three target hull regimes.
    for index in [0,12,21]:
        m,f=target('888',index);key=f'888_{index}_20';cert=json.loads((HERE/f'certificates/{key}.json').read_text())
        a=cpp(m,f,mu=6,counts=(2,0),binary=OLD,nodes=100,seconds=.2,leaves=3)
        for xs in a['leaves']:
            normal=joint_normal(xs,cert)
            b=cpp(m,f,mu=6,counts=(2,0),prefix=normal,symmetry=HERE/f'inputs/{key}.sym',nodes=8,seconds=.3)
            assert b['complete'] and b['leaves']==[normal]
        report['bounded_prefixes'].append(dict(class_index=index,hull=m.hull,stage3=a,covered=len(a['leaves']),scope='Bounded target root sample'))
    for limits in [dict(nodes=0),dict(nodes=1),dict(leaves=0),dict(seconds=0)]:
        r=cpp(m,f,mu=6,**limits);assert r['status']=='LIMIT' and not r['complete']
    report['limits']='PASS';report['seconds']=time.monotonic()-started
    report['hashes']=hashes();atomic_json(HERE/'verification.json',report)
    print(json.dumps(dict(status='PASS',fixtures=len(report['fixtures']),target_corpus=len(report['target_corpus']),
        sampled_target_leaves=sum(r['covered'] for r in report['bounded_prefixes']),positive_controls=report['positive_controls'],seconds=report['seconds'])),flush=True)


if __name__=='__main__':main()
