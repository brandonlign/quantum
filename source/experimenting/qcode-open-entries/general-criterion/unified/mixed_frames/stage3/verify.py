"""Exact small controls, disabled-prune equality and saved-evidence regressions."""
from collections import Counter
import copy
import itertools
import json
from pathlib import Path
import time
from model import ExtensionModel,dual,embed,gram,lift_coordinate,span,weight,symp
from mixed_frames.frames import adapted_frame,xor_table
from mixed_frames.verify import independent_xs,brute_profile
from mixed_frames.stage2.common import BINARY as OLD_BINARY
from mixed_frames.stage2.verify import profile_certificates
from .common import HERE,ROOT,atomic_json,graph_counts,physical_translate_costs,target,budget,sha,input_text
from .checks import cpp
from .run import complete_leaf,check_completion,parse_log
from . import run
import tempfile

def physical_counts(model,frame,columns,w6=6,w7=7):
    converted=frame.model_columns(model,columns)
    assert model.criterion(converted)
    S,_=model.construct(converted)
    assert all(symp(a,b,model.n)==0 for a in S for b in S)
    R={embed(r,0,model.m,model.f) for r in span(frame.R)}
    hist=Counter(weight(w,model.n) for w in span(S) if w not in R)
    return dict(minimum=min(hist),light=sum(n for w,n in hist.items() if w<=w7),n6=hist[w6],n7=hist[w7],histogram={str(w):n for w,n in sorted(hist.items())})

def reject(fn):
    try:fn()
    except (AssertionError,KeyError):return
    raise AssertionError('Malformed evidence accepted')

def main():
    started=time.monotonic();report=dict(status='PASS',equality=[],controls=[],integrity=[])
    # Reuse the independent integer verifier with its sole output redirected.
    profile_certificates.__globals__['HERE']=HERE
    profile_certificates()
    fixtures=[ExtensionModel(6,0,3,4,[15,240],[]),ExtensionModel(7,1,2,4,[240,15],[24,3])]
    for m in fixtures:
        f=adapted_frame(4,m.R)
        for rank,domains in itertools.product([False,True],repeat=2):
            for mu in (0,2,3,4):
                kw=dict(mu=mu,rank=rank,domains=domains)
                a=cpp(m,f,**kw,binary=OLD_BINARY);b=cpp(m,f,**kw)
                assert (a['leaves'],a['nodes'],a['status'])==(b['leaves'],b['nodes'],b['status'])
                report['equality'].append(dict(n=m.n,mu=mu,rank=rank,domains=domains,leaves=b['leaf_count']))
        # All compatible physical graphs; known Y controls need no solver assumptions.
        total=0;positive=0
        for xs in independent_xs(f.p,m.b):
            compatible,_=brute_profile(m,f,xs)
            for ys in compatible:
                columns=list(xs)+list(ys);low=2 if m.hull else 3
                c=physical_counts(m,f,columns,low,low+1)
                assert graph_counts(m,f,columns,low,low+1)==c
                if c['minimum']<low:continue
                counts=(c['n6'],c['n7'])
                r=cpp(m,f,mu=low,counts=counts,w6=low,w7=low+1,prefix=xs)
                assert r['status']=='COMPLETE' and r['leaves']==[list(xs)]
                positive+=1
            total+=len(compatible)
        report['controls'].append(dict(n=m.n,hull=m.hull,physical_graphs=total,positive_cpp=positive))
    # Singular Gamma, sparse physical S, active *proper-depth* mu and count bounds.
    # S has only 16 words: suffix P={Z0,Z1Z2Z3,Z4Z5Z6}; R=Z^4.
    P=[1<<7,14<<7,112<<7]
    m=ExtensionModel(11,7,1,4,[240],dual(P,7));f=adapted_frame(4,m.R)
    assert m.hull==3 and not any(m.Gamma)
    # Relabel the intended coefficient map to the actual independent P basis.
    from model import coordinates
    intended=[1,3,5,0,0,0]
    change=[coordinates(h,P) for h in m.P]
    columns=[sum(((a&c).bit_count()&1)<<i for i,a in enumerate(change)) for c in intended]
    measured=physical_counts(m,f,columns,5,6)
    assert graph_counts(m,f,columns,5,6)==measured and measured['minimum']==5
    for j in (1,2):
        q=1<<f.p
        caps=[max(f.costs[z|(x<<f.p)] for x in range(q) if x&((1<<j)-1)==low for z in range(q)) for low in range(1<<j)]
        active=[a for a,w in enumerate([weight(h,m.f) for h in xor_table(m.P)]) if a and w+min(caps)<5]
        assert active
    r=cpp(m,f,mu=5,counts=(measured['n6'],measured['n7']),w6=5,w7=6,prefix=columns[:3])
    assert r['status']=='COMPLETE' and r['leaves']==[columns[:3]]
    # With insufficient allowance an active partial count test must reject.
    denied=cpp(m,f,mu=5,counts=(0,0),w6=5,w7=6,prefix=columns[:3])
    assert denied['status']=='COMPLETE' and not denied['leaves'] and denied['count_prunes']>0
    # We also require the proper-depth mu test to actively reject a bad X prefix.
    # P ordering may put the singleton at another coefficient; pick it explicitly below.
    singleton=next(a for a,h in enumerate(xor_table(m.P)) if weight(h,7)==1)
    badbeta=next(c for c in range(1,8) if not ((c&singleton).bit_count()&1))
    bad=cpp(m,f,mu=5,prefix=[badbeta],depth=1,rank=False,domains=False)
    assert not bad['leaves'] and bad['mu_prunes']>0
    report['controls'].append(dict(n=11,hull=m.hull,columns=columns,measured=measured,
         active_depths=[1,2],survives=True,denied_count_prunes=denied['count_prunes'],bad_mu_prunes=bad['mu_prunes']))
    # Same construction at the actual weight-6/7 thresholds (16 physical S words).
    P=[3<<10,60<<10,960<<10]
    m=ExtensionModel(14,10,1,4,[240],dual(P,10));f=adapted_frame(4,m.R)
    change=[coordinates(h,P) for h in m.P]
    columns=[sum(((a&c).bit_count()&1)<<i for i,a in enumerate(change)) for c in intended]
    measured=physical_counts(m,f,columns)
    assert measured==graph_counts(m,f,columns) and measured['minimum']==6
    r=cpp(m,f,mu=6,counts=(measured['n6'],measured['n7']),prefix=columns[:3])
    assert r['leaves']==[columns[:3]]
    # Ann(X) knows exact weights of every R translate, not just one minimum.
    small=fixtures[0];sf=adapted_frame(4,small.R)
    exact=cpp(small,sf,mu=3,counts=(10,0),w6=3,w7=4)
    assert exact['status']=='COMPLETE' and exact['exact_prunes']>0 and exact['count_prunes']==0
    report['controls'].append(dict(n=14,hull=m.hull,columns=columns,measured=measured,
          actual_weight6_weight7_thresholds=True,survives=True,exact_translate_prunes=exact['exact_prunes']))
    # Generalized profile completion pipeline, exact equality after Y.
    m=fixtures[0];f=adapted_frame(4,m.R);xs=[3,12]
    base=complete_leaf(m,f,xs,3,None,3,4)
    expected={'fixture':'v1'};base['hashes']=expected
    check_completion(base,m,f,xs,3,None,expected,3,4)
    assert base['accepted']
    for counts in {(c['n6'],c['n7']) for c in base['measurements']}:
        done=complete_leaf(m,f,xs,3,counts,3,4);done['hashes']=expected
        check_completion(done,m,f,xs,3,counts,expected,3,4)
        assert done['accepted']==[s for s,c in zip(base['completion']['solutions'],base['measurements']) if c['minimum']>=3 and (c['n6'],c['n7'])==counts]
    reject(lambda:check_completion(base,m,f,xs,4,None,expected,3,4))
    reject(lambda:check_completion(base,m,f,xs,3,(0,0),expected,3,4))
    corrupt=copy.deepcopy(base);corrupt['accepted']=[]
    reject(lambda:check_completion(corrupt,m,f,xs,3,None,expected,3,4))
    summary=dict(type='summary',status='COMPLETE',complete=True,leaf_count=1)
    cp=dict(enumerator=summary,leaves=1)
    log=json.dumps(dict(type='leaf',X=xs))+'\n'+json.dumps(summary)+'\n'
    assert parse_log(log,cp,2)==[xs]
    for text in ['',json.dumps(dict(type='leaf',X=xs)),log+json.dumps(summary)+'\n']:
        reject(lambda text=text:parse_log(text,cp,2))
    corrupt=copy.deepcopy(cp);corrupt['enumerator']['leaf_count']=0
    reject(lambda:parse_log(log,corrupt,2))
    # Actual saved-partition verifier rejects the two aggregate audit attacks.
    row=dict(frame='fixture',class_index=0,profile=[0,0],budget=None,d=m.d,mu_out=3,input_sha256='fixture')
    evidence=complete_leaf(m,f,xs,3);evidence['hashes']=expected
    assert evidence['accepted']
    old_target=run.target
    run.target=lambda name,index:(m,f)
    try:
        with tempfile.TemporaryDirectory(dir=HERE,prefix='integrity_') as folder_name:
            folder=Path(folder_name);(folder/'completions').mkdir()
            name=run.key(xs);atomic_json(folder/'completions'/name,evidence)
            (folder/'enumerator.jsonl').write_text(log)
            cp=dict(**run.context(row,xs[0],expected),status='COMPLETE',complete=True,excluded=False,
                    enumerator=summary,leaves=1,rejected_leaves=0,surviving_solutions=len(evidence['accepted']),
                    log_sha256=sha(folder/'enumerator.jsonl'),completion_sha256={name:sha(folder/'completions'/name)})
            atomic_json(folder/'checkpoint.json',cp)
            assert run.verify_partition(folder,row,xs[0],expected)['surviving_solutions']>0
            attack=copy.deepcopy(cp);attack['surviving_solutions']=0;attack['excluded']=True
            atomic_json(folder/'checkpoint.json',attack)
            reject(lambda:run.verify_partition(folder,row,xs[0],expected))
            (folder/'enumerator.jsonl').write_text('');attack=copy.deepcopy(cp)
            attack['log_sha256']=sha(folder/'enumerator.jsonl');atomic_json(folder/'checkpoint.json',attack)
            reject(lambda:run.verify_partition(folder,row,xs[0],expected))
        # Run the real worker/ack pipeline with a poisoned mu=5 cache, then
        # request mu=4. The physical n=6 regression must recover its solutions.
        old_here=run.HERE
        try:
            with tempfile.TemporaryDirectory(dir=HERE,prefix='worker_') as folder_name:
                dest=Path(folder_name);run.HERE=dest
                inp=dest/'fixture.txt';inp.write_text(input_text(m,f,4))
                worker_row=dict(row,key='fixture',mu_out=4,input='fixture.txt',input_sha256=sha(inp))
                cache=dest/'runs'/'fixture'/str(xs[0])/'completions';cache.mkdir(parents=True)
                poisoned=complete_leaf(m,f,xs,5);poisoned['hashes']=expected
                assert not poisoned['accepted']
                atomic_json(cache/run.key(xs),poisoned)
                answer=run.partition((worker_row,xs[0],expected))
                assert answer['status']=='COMPLETE' and not answer['excluded'] and answer['surviving_solutions']>0,answer
                restored=json.loads((cache/run.key(xs)).read_text())
                assert restored['context']['mu_out']==4 and len(restored['accepted'])==2
                report['worker_cache_regression']=dict(status=answer['status'],surviving_solutions=answer['surviving_solutions'],cached_leaf_recovered_solutions=2)
        finally:run.HERE=old_here
    finally:run.target=old_target
    report['integrity']=['missing/multiple/mismatching terminal summary rejected','cached mu/profile mismatch rejected','accepted list and partition survival totals recomputed exactly']
    # Bounded target comparisons in every hull regime and both frames.
    for name,index in [('688',7),('688',0),('688',12),('888',7),('888',0),('888',12)]:
        m,f=target(name,index)
        for mu in (0,6,7,8):
            options=dict(mu=mu,nodes=12,leaves=3)
            a=cpp(m,f,**options,binary=OLD_BINARY);b=cpp(m,f,**options)
            assert (a['leaves'],a['nodes'],a['status'])==(b['leaves'],b['nodes'],b['status'])
            report['equality'].append(dict(frame=name,index=index,hull=m.hull,mu=mu,nodes=b['nodes'],leaves=b['leaf_count']))
    # Nineteen actual mu6 necessary Z frontiers from the old 888 screen.
    # Their exact known pure-Z graph words are checked physically, not with caps.
    from .frontier_pilot import legacy_frontiers
    corpus=[]
    for item in legacy_frontiers():
        m,f=target('888',item['class_index']);xs=item['X']
        old=cpp(m,f,mu=6,prefix=xs,binary=OLD_BINARY)
        unlimited=cpp(m,f,mu=6,prefix=xs)
        assert old['leaves']==unlimited['leaves']==[xs]
        words=xor_table(m.P);oldZ=[((1<<i)|(1<<7))<<8 for i in range(6)]
        exact_hist=Counter()
        for z,a in enumerate(xor_table(item['legacy_Z'])):
            if not z:continue
            for r in xor_table(f.R):
                exact_hist[weight(lift_coordinate(z,oldZ)^r,8)+weight(words[a],9)]+=1
        # At full depth, Ann(X) labels these same 63 physical suffixes.
        r=cpp(m,f,mu=6,counts=(2,0),prefix=xs)
        necessary=exact_hist[6]<=2 and exact_hist[7]==0
        assert bool(r['leaves'])==necessary
        corpus.append(dict(**item,known_n6=exact_hist[6],known_n7=exact_hist[7],
                           baseline_leaf=True,survives_profile20=bool(r['leaves']),
                           count_prunes=r['count_prunes'],exact_prunes=r['exact_prunes'],proper_exact_prunes=r['proper_exact_prunes']))
    assert any(r['proper_exact_prunes'] for r in corpus)
    report['target_corpus']=corpus
    for kw in [dict(nodes=0),dict(nodes=1),dict(leaves=0),dict(leaves=1),dict(seconds=0)]:
        r=cpp(fixtures[0],adapted_frame(4,fixtures[0].R),**kw)
        assert r['status']=='LIMIT' and not r['complete']
    report['limits']='PASS'
    report['seconds']=time.monotonic()-started
    report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'enumerate_z.cpp',HERE/'enumerate_z',HERE/'verify.py',OLD_BINARY]}
    atomic_json(HERE/'verification.json',report)
    print(json.dumps(dict(status='PASS',equality=len(report['equality']),controls=report['controls'],seconds=report['seconds'])),flush=True)

if __name__=='__main__':main()
