"""Serial exhaustive pipeline, explicitly launched with --exhaustive.

Independent saved-evidence verification is also used before COMPLETE reuse.
All I/O is under stage3; legacy stage-2 evidence is never reused as a cache.
"""
import argparse
import copy
import hashlib
import json
import os
import subprocess
import time
import traceback
from mixed_frames.fixed_z import solve
from .common import (HERE,BINARY,atomic_json,budget,case_input,graph_counts,
                     hashes,sha,target)

def key(xs):return hashlib.sha256(json.dumps(xs).encode()).hexdigest()+'.json'

def completion_context(model,frame,mu,counts,w6=6,w7=7):
    return dict(n=model.n,k=model.k,d=model.d,m=model.m,R=frame.R,E=frame.E,
                C=model.C,P=model.P,Gamma=model.Gamma,hull=model.hull,
                mu_out=mu,budget=list(counts) if counts is not None else None,w6=w6,w7=w7)

def accepted(count,mu,counts):
    return count['minimum']>=mu and (counts is None or (count['n6'],count['n7'])==tuple(counts))

def complete_leaf(model,frame,xs,mu,counts=None,w6=6,w7=7):
    result=solve(model,frame,xs,find_all=True)
    assert result['status']=='COMPLETE' and result['complete']
    measurements=[graph_counts(model,frame,s['adapted_columns'],w6,w7) for s in result['solutions']]
    keep=[copy.deepcopy(s) for s,c in zip(result['solutions'],measurements) if accepted(c,mu,counts)]
    return dict(X=list(xs),status='COMPLETE',complete=True,
                context=completion_context(model,frame,mu,counts,w6,w7),
                raw_solutions=len(result['solutions']),measurements=measurements,
                accepted=keep,completion=result)

def check_completion(c,model,frame,xs,mu,counts,expected,w6=6,w7=7):
    assert c['hashes']==expected and c['X']==list(xs)
    assert c['context']==completion_context(model,frame,mu,counts,w6,w7), 'Cached completion context changed'
    assert c['status']=='COMPLETE' and c['complete']
    result=c['completion']
    assert result['status']=='COMPLETE' and result['complete'] and result['X']==list(xs)
    assert result['limits']==dict(nodes=None,seconds=None,solutions=None)
    ctx=result['context']
    for name,val in dict(n=model.n,k=model.k,d=model.d,m=model.m,R=frame.R,E=frame.E,C=model.C,P=model.P,hull=model.hull,frame=frame.name).items():
        assert ctx[name]==val
    assert c['raw_solutions']==len(result['solutions'])
    columns=[tuple(s['adapted_columns']) for s in result['solutions']]
    assert len(set(columns))==len(columns)
    assert all(list(s[:frame.p])==list(xs) for s in columns)
    assert all(model.criterion(frame.model_columns(model,list(s))) for s in columns), 'Invalid cached physical map'
    measurements=[graph_counts(model,frame,s['adapted_columns'],w6,w7) for s in result['solutions']]
    assert c['measurements']==measurements
    assert c['accepted']==[s for s,m in zip(result['solutions'],measurements) if accepted(m,mu,counts)]
    return len(c['accepted'])

def parse_log(text,checkpoint,p):
    rows=[json.loads(line) for line in text.splitlines()]
    assert rows and rows[-1]['type']=='summary', 'Missing terminal summary'
    assert sum(r['type']=='summary' for r in rows)==1, 'Multiple summaries'
    summary=rows[-1]
    assert summary==checkpoint['enumerator'], 'Terminal/checkpoint summary mismatch'
    assert summary['status']=='COMPLETE' and summary['complete']
    assert all(r['type']=='leaf' and len(r['X'])==p for r in rows[:-1])
    xs=[r['X'] for r in rows[:-1]]
    assert len(xs)==summary['leaf_count']==checkpoint['leaves']
    assert len({tuple(x) for x in xs})==len(xs)
    return xs

def context(row,first,expected):
    return dict(frame=row['frame'],class_index=row['class_index'],profile=row['profile'],
                budget=row['budget'],d=row['d'],mu_out=row['mu_out'],first_image=first,
                input_sha256=row['input_sha256'],hashes=expected)

def verify_partition(folder,row,first,expected):
    cp=json.loads((folder/'checkpoint.json').read_text())
    assert all(cp[k]==v for k,v in context(row,first,expected).items())
    assert cp['status']=='COMPLETE' and cp['complete']
    model,frame=target(row['frame'],row['class_index'])
    assert sha(folder/'enumerator.jsonl')==cp['log_sha256']
    xs=parse_log((folder/'enumerator.jsonl').read_text(),cp,frame.p)
    assert all(x[0]==first for x in xs)
    assert set(cp['completion_sha256'])=={key(x) for x in xs}
    surviving=rejected=0
    for x in xs:
        name=key(x);path=folder/'completions'/name
        assert sha(path)==cp['completion_sha256'][name]
        n=check_completion(json.loads(path.read_text()),model,frame,x,row['mu_out'],row['budget'],expected)
        surviving+=n;rejected+=n==0
    assert cp['surviving_solutions']==surviving and cp['rejected_leaves']==rejected
    assert cp['excluded']==(surviving==0)
    return cp

def partition(task):
    row,first,expected=task
    folder=HERE/'runs'/row['key']/str(first);folder.mkdir(parents=True,exist_ok=True)
    path=folder/'checkpoint.json'
    if path.exists():
        saved=json.loads(path.read_text())
        assert all(saved[k]==v for k,v in context(row,first,expected).items())
        if saved['status']=='COMPLETE':return verify_partition(folder,row,first,expected)
    assert sha(HERE/row['input'])==row['input_sha256']
    model,frame=target(row['frame'],row['class_index'])
    count=rejected=surviving=0;done_files={};start=time.monotonic();proc=None
    completion_dir=folder/'completions';completion_dir.mkdir(exist_ok=True)
    def save(status,**extra):
        v=dict(**context(row,first,expected),status=status,complete=status=='COMPLETE',
               excluded=status=='COMPLETE' and surviving==0,leaves=count,rejected_leaves=rejected,
               surviving_solutions=surviving,seconds=time.monotonic()-start,**extra)
        atomic_json(path,v);return v
    save('RUNNING')
    try:
        with (folder/'stderr.log').open('w') as err,(folder/'enumerator.jsonl').open('w') as log:
            proc=subprocess.Popen([str(BINARY),str(HERE/row['input']),'--prefix',str(first),'--ack','1'],
                                  stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
            summary=None
            for line in proc.stdout:
                log.write(line);log.flush();obj=json.loads(line)
                assert summary is None, 'Output after terminal summary'
                if obj['type']=='summary':summary=obj;continue
                assert obj['type']=='leaf' and len(obj['X'])==frame.p
                xs=obj['X'];name=key(xs);assert name not in done_files
                file=completion_dir/name
                if file.exists():
                    done=json.loads(file.read_text())
                    # Reject a stale context; never count it as an exclusion.
                    try:check_completion(done,model,frame,xs,row['mu_out'],row['budget'],expected)
                    except (AssertionError,KeyError):done=None
                else:done=None
                if done is None:
                    done=complete_leaf(model,frame,xs,row['mu_out'],row['budget'])
                    done['hashes']=expected;atomic_json(file,done)
                n=check_completion(done,model,frame,xs,row['mu_out'],row['budget'],expected)
                count+=1;surviving+=n;rejected+=n==0;done_files[name]=sha(file)
                save('RUNNING');proc.stdin.write('OK\n');proc.stdin.flush()
            assert proc.wait()==0 and summary is not None
            assert summary['leaf_count']==count
            assert summary['complete']==(summary['status']=='COMPLETE')
            status='COMPLETE' if summary['complete'] else 'LIMIT'
        saved=save(status,enumerator=summary,log_sha256=sha(folder/'enumerator.jsonl'),completion_sha256=done_files)
        if status=='COMPLETE':return verify_partition(folder,row,first,expected)
        return saved
    except BaseException:
        if proc is not None and proc.poll() is None:proc.terminate();proc.wait()
        return save('ERROR',error=traceback.format_exc())

def read_manifest():
    manifest=json.loads((HERE/'manifest.json').read_text())
    assert hashes()==manifest['hashes'],'Implementation changed; prepare again before a new census'
    for name,digest in manifest['table_sha256'].items():assert sha(HERE/name)==digest
    records=json.loads((HERE/'capacity_certificate.json').read_text())['records']
    from .common import PROFILES,prefilter,profile_mu
    assert len(records)==112 and {(r['frame'],r['class_index']) for r in records}=={(name,i) for name in PROFILES for i in range(56)}
    pairs={(r['frame'],r['class_index']) for r in records if r['verdict']=='CANDIDATE'}
    required={(name,i,tuple(p)) for name,i in pairs for p in PROFILES[name]}
    assert {(r['frame'],r['class_index'],tuple(r['profile'])) for r in manifest['classes']}==required
    assert len(manifest['classes'])==len(required)
    for rec in records:
        m,f=target(rec['frame'],rec['class_index'])
        assert rec==dict(frame=rec['frame'],class_index=rec['class_index'],**prefilter(m,f))
    for row in manifest['classes']:
        m,f=target(row['frame'],row['class_index'])
        assert row['budget']==list(budget(f,row['profile'])) and row['d']==m.d==6
        assert row['mu_out']==profile_mu(f,row['profile'])
        for k,v in dict(R=f.R,E=f.E,C=m.C,P=m.P,Gamma=m.Gamma,hull=m.hull).items():assert row[k]==v
        assert row['first_images']==[s for s in range(1<<m.b) if s and m.distances[s]+f.costs[1]>=6]
        assert (HERE/row['input']).read_text()==case_input(m,f,row['profile'])
        assert sha(HERE/row['input'])==row['input_sha256']
    return manifest

def aggregate():
    manifest=read_manifest();classes=[]
    for row in manifest['classes']:
        parts=[]
        for first in row['first_images']:
            folder=HERE/'runs'/row['key']/str(first);path=folder/'checkpoint.json'
            if not path.exists():r=dict(status='PENDING',complete=False,excluded=False)
            else:
                r=json.loads(path.read_text())
                assert all(r[k]==v for k,v in context(row,first,manifest['hashes']).items())
                assert r['complete']==(r['status']=='COMPLETE')
                if r['complete']:r=verify_partition(folder,row,first,manifest['hashes'])
                else:assert not r['excluded']
            parts.append(r)
        complete=all(r['complete'] for r in parts)
        classes.append(dict(key=row['key'],frame=row['frame'],class_index=row['class_index'],profile=row['profile'],
                            complete=complete,excluded=complete and all(r['excluded'] for r in parts),
                            complete_partitions=sum(r['complete'] for r in parts),partitions=len(parts)))
    result=dict(status='COMPLETE' if all(c['complete'] for c in classes) else 'INCOMPLETE',verified=True,
                classes=classes,manifest_sha256=sha(HERE/'manifest.json'))
    atomic_json(HERE/'aggregate.json',result);return result

def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exhaustive',action='store_true')
    parser.add_argument('--aggregate',action='store_true')
    parser.add_argument('--frames',nargs='+',choices=['688','888'],default=['688','888'])
    parser.add_argument('--classes',nargs='+',type=int)
    args=parser.parse_args()
    if args.aggregate:print(json.dumps(aggregate()));return
    if not args.exhaustive:parser.error('Explicit --exhaustive is required; pilots use pilot.py')
    manifest=read_manifest()
    for row in manifest['classes']:
        if row['frame'] not in args.frames or args.classes is not None and row['class_index'] not in args.classes:continue
        for first in row['first_images']:
            r=partition((row,first,manifest['hashes']))
            print(json.dumps({k:r[k] for k in ('frame','class_index','profile','first_image','status','leaves','seconds')}),flush=True)
            if r['status']=='ERROR':raise RuntimeError(r['error'])
    aggregate()

if __name__=='__main__':main()
