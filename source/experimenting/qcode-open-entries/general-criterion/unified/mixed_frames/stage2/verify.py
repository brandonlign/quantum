"""Read-only stage-1 replay plus stage-2 exact leaf, graph, pipeline tests.

Only artifacts under stage2/ are written. Run with -B, assertions enabled.
"""
import contextlib
import io
import itertools
import json
import runpy
import subprocess
import time
from pathlib import Path

from model import ExtensionModel, embed, lift_coordinate, reduced_basis, span, weight
from mixed_frames.frames import Frame, adapted_frame, mixed_frame, xor_table
from mixed_frames.enumerate_z import enumerate_z
from mixed_frames.fixed_z import solve, validate_solution
from mixed_frames.verify import positive_controls, brute_profile, independent_xs, physical_weight
from mixed_frames.suffix_tables import replay
from .common import HERE, ROOT, BINARY, atomic_json, graph_minimum, graph_z_ok, graph_partial_ok, sha, target, write_input
from .reference import enumerate_reference
from .run import complete_leaf


def cpp(model,frame,*,mu=0,prefix=(),depth=None,node_limit=None,leaf_limit=None,rank=True,domains=True,seconds=None):
    path=HERE/'validation_input.txt'
    write_input(path,model,frame,mu,rank,domains)
    cmd=[str(BINARY),str(path)]
    for option,value in [('--prefix',','.join(map(str,prefix)) if prefix else None),('--depth',depth),
                         ('--nodes',node_limit),('--leaves',leaf_limit),('--seconds',seconds)]:
        if value is not None: cmd += [option,str(value)]
    proc=subprocess.run(cmd,text=True,capture_output=True,check=True)
    rows=[json.loads(line) for line in proc.stdout.splitlines()]
    assert rows[-1]['type']=='summary'
    return dict(**rows[-1],leaves=[r['X'] for r in rows[:-1]])


def compare(model,frame,**kw):
    refkw=dict(kw)
    if 'rank' in refkw:refkw['use_rank']=refkw.pop('rank')
    if 'domains' in refkw:refkw['use_domains']=refkw.pop('domains')
    a=enumerate_reference(model,frame,**refkw)
    c=cpp(model,frame,**kw)
    assert a['leaves']==c['leaves'], (frame.name,kw,a,c)
    assert a['status']==c['status'] and a['nodes']==c['nodes'], (kw,a,c)
    assert len({tuple(x) for x in c['leaves']})==len(c['leaves'])
    return dict(frame=frame.name,parameters=kw,status=c['status'],nodes=c['nodes'],leaves=c['leaf_count'],
                mu_prunes=c['mu_prunes'])


def profile_certificates():
    checker=runpy.run_path(str(ROOT.parents[1]/'pure/check_general.py'))
    bundle=json.loads((ROOT/'weight8_split/frame_count_certificates.json').read_text())
    class MemoryPath:
        def __init__(self,name):self.name=name
        def read_text(self):return json.dumps(bundle[self.name])
    checker['check'].__globals__['Path']=MemoryPath
    log=io.StringIO()
    with contextlib.redirect_stdout(log):
        for name in bundle:checker['check'](name)
    assert len(bundle)==13
    # Independent ordinary-enumerator row reconstruction; no scipy or LP.
    kraw=checker['kraw'];M=16384
    residual=[]
    for j in range(18):
        v=-1835008*int(j==0)+116
        for w,mul in enumerate([60,25,11,3,1],1):v+=mul*(kraw(17,w,j)-M*int(j==w))
        for w,mul in enumerate([-851968,-278528,-507904,-147456,-180224],1):v+=mul*int(j==w)
        for w,mul in enumerate([88,32,18,4,2]):v-=mul*((-1)**j)*kraw(17,w,j)
        residual.append(v-32768*int(j in (6,7)))
    assert residual==[0]*18 and -1835008+116*M==65536
    atomic_json(HERE/'profile_verification.json',dict(status='PASS',certificates=13,
                log=log.getvalue(),global_weight_sum_identity_residual=residual,rhs=65536,
                bundle_sha256=sha(ROOT/'weight8_split/frame_count_certificates.json'),
                checker_sha256=sha(ROOT.parents[1]/'pure/check_general.py')))


def physical_mu(model,frame,columns):
    S,_=model.construct(frame.model_columns(model,columns))
    R={embed(r,0,model.m,model.f) for r in span(frame.R)}
    return min(weight(w,model.n) for w in span(S) if w not in R)


def main():
    started=time.monotonic();report=dict(status='PASS',comparisons=[],controls=[],pipeline_cases=0)
    profile_certificates()
    print('PASS 13 frame certificates and exact A6+A7 identity',flush=True)
    # The existing replay writes nothing; includes all physical suffix tables.
    replay()
    print('PASS stage-1 capacity/prefix/suffix physical replay',flush=True)
    controls=positive_controls()
    assert len(controls)==7
    report['stage1_positive_controls']=[c['label'] for c in controls]
    print('PASS stage-1 positive physical controls',flush=True)
    fixtures=[(ExtensionModel(6,0,3,4,[15,240],[]),None),
              (ExtensionModel(7,1,3,4,[240,15],[24,3]),None)]
    for model,_ in fixtures:
        frame=adapted_frame(4,model.R)
        for rank,domains in itertools.product([False,True],repeat=2):
            ref=enumerate_z(model,frame,use_rank=rank,use_domains=domains)
            row=compare(model,frame,rank=rank,domains=domains)
            got=cpp(model,frame,rank=rank,domains=domains)
            assert got['leaves']==ref['leaves'] and got['nodes']==ref['nodes']
            report['comparisons'].append(row)
        for mu in (3,4,5):
            report['comparisons'].append(compare(model,frame,mu=mu))
        # Every compatible small graph: exact physical S minimum, and no
        # false partial or pure-Z prune when its measured mu is requested.
        physical_graphs=0
        for xs in independent_xs(frame.p,model.b):
            profile,_=brute_profile(model,frame,xs)
            for ys in profile:
                columns=list(xs)+list(ys)
                measured=physical_mu(model,frame,columns)
                assert graph_minimum(model,frame,columns)==measured
                assert graph_z_ok(model,frame,list(xs),measured)
                assert all(graph_partial_ok(model,frame,list(xs[:j]),measured) for j in range(frame.p+1))
                physical_graphs+=1
        report.setdefault('physical_graph_checks',[]).append(dict(n=model.n,graphs=physical_graphs))
        # Pipeline exact accepted solution sets at nonempty and rejecting mu.
        for xs in independent_xs(frame.p,model.b)[:8]:
            baseline=solve(model,frame,list(xs),find_all=True)
            for mu in (3,4,5):
                done=complete_leaf(model,frame,list(xs),mu)
                expected={tuple(s['adapted_columns']) for s in baseline['solutions']
                          if physical_mu(model,frame,s['adapted_columns'])>=mu}
                assert {tuple(s['adapted_columns']) for s in done['accepted']}==expected
                report['pipeline_cases']+=1
    print('PASS small exact Z leaf lists, independent graph weights, exhaustive completion pipeline',flush=True)
    # 35 length-two subspaces exercise all hull regimes in a complete census.
    subspaces={tuple(reduced_basis([a,b])) for a in range(1,16) for b in range(a+1,16)
               if len(reduced_basis([a,b]))==2}
    for C in sorted(subspaces):
        model=ExtensionModel(4,1,2,2,[12],list(C));frame=adapted_frame(2,model.R)
        for d in (2,3):
            model.d=d
            ref=enumerate_z(model,frame)
            got=cpp(model,frame)
            assert got['status']=='COMPLETE' and got['leaves']==ref['leaves']
            report['comparisons'].append(dict(n=4,C=list(C),d=d,leaves=got['leaf_count']))
    # Real six-column mixed controls: mu is measured in S, not D.
    for name in ('668','778','688'):
        frame=mixed_frame(name);model=ExtensionModel(14,0,2,8,frame.R,[])
        columns=[model.syndrome(1<<(6+i)) for i in range(6)]+[model.syndrome(1<<i) for i in range(6)]
        validate_solution(model,frame,columns)
        mu=physical_mu(model,frame,columns)
        assert graph_minimum(model,frame,columns)==mu
        result=cpp(model,frame,mu=mu,prefix=columns[:6])
        assert result['status']=='COMPLETE' and result['leaves']==[columns[:6]]
        # This known Y passes the exact final predicate. Exhaustive pipeline
        # comparisons above use the small nonempty fixtures; enumerating all
        # 2^21 d=2 six-variable completions is unnecessary for this control.
        assert graph_minimum(model,frame,columns)>=mu
        report['controls'].append(dict(frame=name,n=14,d=2,mu_out=mu,known_Y_survives=True))
    print('PASS mixed-frame measured-mu positive controls',flush=True)
    # Bounded target prefixes: actual labels, hull 0/2/6, frontiers and full
    # traversal; compare ordered leaves, status and exact visited-node counts.
    for name,index in [('668',0),('668',7),('778',7),('778',12),('688',7)]:
        model,frame=target(name,index)
        first=next(s for s,d in enumerate(model.distances) if d>=5)
        for mu in (0,8 if name!='688' else 6):
            for options in [dict(prefix=[first],depth=2,leaf_limit=12),
                            dict(prefix=[first],node_limit=25,leaf_limit=8)]:
                row=compare(model,frame,mu=mu,**options);row['class_index']=index
                report['comparisons'].append(row)
        # Unpartitioned original reference versus C++, mu disabled.
        a=enumerate_z(model,frame,node_limit=12,leaf_limit=3)
        c=cpp(model,frame,node_limit=12,leaf_limit=3)
        assert a['leaves']==c['leaves'] and a['status']==c['status'] and a['nodes']==c['nodes']
    print('PASS bounded target prefixes and original reference traversal',flush=True)
    model,frame=fixtures[0][0],adapted_frame(4,fixtures[0][0].R)
    for options in [dict(node_limit=0),dict(node_limit=1),dict(leaf_limit=0),dict(leaf_limit=1),dict(seconds=0)]:
        r=cpp(model,frame,**options)
        assert r['status']=='LIMIT' and not r['complete']
    report['limits']='PASS nodes/time/zero and nonzero leaf caps; none is COMPLETE'
    # Exercise the real leaf/ack protocol with exhaustive completion between
    # emission and acknowledgement, including LIMIT after an emitted leaf.
    write_input(HERE/'validation_input.txt',model,frame)
    proc=subprocess.Popen([str(BINARY),str(HERE/'validation_input.txt'),'--leaves','1','--ack','1'],
                          stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    leaf=json.loads(proc.stdout.readline())
    assert leaf['type']=='leaf'
    done=complete_leaf(model,frame,leaf['X'],3)
    assert done['complete']
    proc.stdin.write('OK\n');proc.stdin.flush()
    summary=json.loads(proc.stdout.readline())
    assert proc.wait()==0 and summary['status']=='LIMIT' and not summary['complete']
    report['ack_pipeline']='PASS exhaustive leaf completion before OK; leaf cap remains LIMIT'
    report['seconds']=time.monotonic()-started
    report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'enumerate_z.cpp',BINARY,
                        HERE/'common.py',HERE/'reference.py',HERE/'verify.py',HERE/'run.py']}
    atomic_json(HERE/'verification.json',report)
    print(json.dumps(dict(status='PASS',comparisons=len(report['comparisons']),seconds=report['seconds'])),flush=True)


if __name__=='__main__':main()
