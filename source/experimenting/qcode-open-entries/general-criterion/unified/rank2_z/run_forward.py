"""Bounded all15 Z+early-X0 runs; completion is recorded explicitly."""
from pathlib import Path
import json,subprocess,time,hashlib,argparse
BASE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=2)
    ap.add_argument('--nodes',type=int,default=100000000);ap.add_argument('--leaf-limit',type=int,default=10)
    ap.add_argument('--classes');ap.add_argument('--name',default='forward_batch');a=ap.parse_args()
    assert json.loads((BASE/'forward_verification.json').read_text())['status']=='PASS'
    inputs=json.loads((BASE.parent/'rank2_z_inputs.json').read_text())['classes'];reports=[];start=time.monotonic()
    chosen=None if a.classes is None else set(map(int,a.classes.split(',')))
    for item in inputs:
        idx=item['class_index']
        if chosen is not None and idx not in chosen:continue
        data=json.loads(Path(item['table']).read_text());m=data['m'];rows=[f'{len(data["flags"])} {m-2} {a.nodes} {int(1000*a.seconds)} {a.leaf_limit}',' '.join(map(str,data['needs']))]
        for flag,form,delta,lam in zip(data['flags'],data['form_rows'],data['coset_distances'],data['suffix_weights']):
            mask=sum(1<<c for c in range(m+1) if delta>=6-c and lam>=8-c)
            assert bool(flag&1)==bool(mask&4) or flag==0
            rows.append(f'{flag} {form} {mask}')
        path=BASE/f'{a.name}_{idx}.txt';path.write_text('\n'.join(rows)+'\n')
        run=subprocess.run([str(BASE/'enumerator_unary'),str(path)],capture_output=True,text=True)
        log=BASE/f'{a.name}_{idx}.jsonl';log.write_text(run.stdout)
        if run.returncode:r=dict(class_index=idx,status='ERROR',Z_excluded=False,error=run.stderr[-2000:])
        else:
            end=json.loads(run.stdout.splitlines()[-1]);r=dict(class_index=idx,**end,log=str(log),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        reports.append(r);print(json.dumps(r),flush=True)
        (BASE/f'{a.name}_checkpoint.json').write_text(json.dumps(dict(reports=reports))+'\n')
    result=dict(status='BATCH_FINISHED',reports=reports,complete_exclusions=sum(r.get('Z_excluded',False) for r in reports),
        limited=sum(r['status']=='LIMIT' for r in reports),surviving_Z_X0_leaves=sum(r.get('leaves',0) for r in reports),
        seconds=time.monotonic()-start,node_limit=a.nodes,time_limit_seconds=a.seconds,leaf_limit=a.leaf_limit,
        scope='Pure-Z distance/isotropy plus necessary X0-image forward domains. '
              'Zero surviving leaves after COMPLETE excludes C; LIMIT is not an exclusion, '
              'surviving leaves need a full X solver.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE/'enumerator_unary.cpp',BASE/'enumerator_unary']})
    (BASE/f'{a.name}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['reports','sha256']}))
if __name__=='__main__':main()
