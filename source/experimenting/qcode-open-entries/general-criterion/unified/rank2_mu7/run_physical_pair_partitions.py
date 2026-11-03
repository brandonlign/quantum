"""Resume verified first-image partitions for the recentered minimum-seven problem."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'opus_parallel/rank_two_extension'))
from rank2_solver import solve

HERE=Path(__file__).resolve().parent
BASE=HERE.parent


def aggregate(records,required):
    if any(r['status']=='FOUND' for r in records):return 'FOUND',[]
    covered=set()
    for r in records:
        if r['status']=='PARTITION_EXCLUDED':covered.update(r['first'])
    return ('EXCLUDED_MIN7_CANDIDATE' if covered==set(required) else 'PARTIAL'),sorted(covered)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--class',dest='index',type=int,required=True)
    ap.add_argument('--seconds',type=float,default=30);ap.add_argument('--chunk',type=int,default=8)
    ap.add_argument('--pair-depth',type=int,default=4);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();assert 0<=a.pair_depth<6 and a.seconds>0 and a.chunk>0 and a.workers in [1,2]
    verification=json.loads((HERE/'physical_pair_verification.json').read_text());assert verification['status']=='PASS'
    for name,digest in verification['sha256'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    source=BASE/f'rank2_mu7_orbits/class_{a.index}.txt';rows=source.read_text().splitlines()
    N,p,*_=map(int,rows[0].split());assert N==4096 and p==6
    required=list(map(int,rows[N+3].split()))
    assert int(rows[N+2])==len(required) and len(required)==len(set(required))
    assert required==sorted(required)
    table=BASE/f'rank2_z/class_{a.index}.json';data=json.loads(table.read_text())
    paths=[source,table,Path(__file__),HERE/'enumerator_physical_pairs.cpp',HERE/'enumerator_physical_pairs',
           HERE/'physical_pair_verification.json',BASE/'rank2_z/aut_orbits.json',
           BASE/'opus_parallel/rank_two_extension/rank2_solver.py']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    records=[]
    if a.output.exists():
        old=json.loads(a.output.read_text())
        assert old['sha256']==hashes and old['class_index']==a.index and old['required_first']==required and old['pair_depth']==a.pair_depth
        records=old['records']
    def save():
        status,covered=aggregate(records,required)
        result=dict(status=status,class_index=a.index,minimum_stabilizer=7,pair_depth=a.pair_depth,required_first=required,
                    covered_first=covered,records=records,sha256=hashes,
                    scope='One nondegenerate recentered min7 candidate only; full target and other candidates remain separate')
        temp=a.output.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(a.output)
    status,covered=aggregate(records,required);save()
    if status in ['FOUND','EXCLUDED_MIN7_CANDIDATE']:print(status);return
    remaining=[r for r in required if r not in covered]
    groups=[remaining[i:i+a.chunk] for i in range(0,len(remaining),a.chunk)]
    artifacts=HERE/f'partitions_{a.index}';artifacts.mkdir(exist_ok=True)
    def run(first):
        key=hashlib.sha256(json.dumps(first).encode()).hexdigest()[:16]
        path=artifacts/f'{key}.txt';log=artifacts/f'{key}.jsonl'
        payload=[f'4096 6 0 {int(a.seconds*1000)} 0']+rows[1:N+2]+[str(len(first)),' '.join(map(str,first))]
        path.write_text('\n'.join(payload)+'\n');verdicts=[];stats=None;found=None
        with log.open('w') as output:
            proc=subprocess.Popen([str(HERE/'enumerator_physical_pairs'),str(path),str(a.pair_depth)],stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,text=True)
            try:
                for line in proc.stdout:
                    output.write(line);output.flush();event=json.loads(line)
                    if event['status']=='LEAF':
                        answer=solve(data['P_generators'],data['C_generators'],8,6,7,event['Z'],node_limit=1000000)
                        verdicts.append(dict(Z=event['Z'],answer=answer))
                        if answer['status']=='FOUND':found=verdicts[-1];proc.terminate();break
                    else:stats=event
                proc.wait();error=proc.stderr.read()
                if found:status='FOUND'
                else:
                    assert proc.returncode==0,error
                    assert stats and stats['status'] in ['COMPLETE','LIMIT']
                    assert stats['leaves']==len(verdicts)
                    status='PARTITION_EXCLUDED' if stats['status']=='COMPLETE' and all(
                        v['answer']['status']=='EXCLUDED' for v in verdicts) else 'LIMIT'
            finally:
                if proc.poll() is None:proc.terminate();proc.wait()
        return dict(status=status,first=first,z=stats,verdicts=verdicts,found=found,
                    input=str(path),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    log=str(log),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
    print(json.dumps(dict(class_index=a.index,pending_first=len(remaining),groups=len(groups))),flush=True)
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(run,group):group for group in groups}
        for future in as_completed(futures):
            r=future.result();records.append(r);save();status,covered=aggregate(records,required)
            print(json.dumps(dict(first=futures[future],status=r['status'],covered=len(covered),total=len(required),
                                  nodes=None if r['z'] is None else r['z']['nodes'])),flush=True)
    print(json.dumps(dict(status=aggregate(records,required)[0],covered=len(aggregate(records,required)[1]))),flush=True)


if __name__=='__main__':main()
