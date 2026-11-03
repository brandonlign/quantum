"""Bounded stage-Z batch: LIMIT and surviving leaves are never exclusions."""
from pathlib import Path
import json,subprocess,hashlib,time
BASE=Path(__file__).resolve().parent

def main():
    checked=json.loads((BASE/'verification.json').read_text());assert checked['status']=='PASS'
    inputs=json.loads((BASE.parent/'rank2_z_inputs.json').read_text())['classes'];reports=[];start=time.monotonic()
    for item in inputs:
        source=Path(item['input']);run=subprocess.run([str(BASE/'enumerator'),str(source)],capture_output=True,text=True)
        log=BASE/f'class_{item["class_index"]}_z.jsonl';log.write_text(run.stdout)
        if run.returncode:
            report=dict(class_index=item['class_index'],status='ERROR',excluded=False,error=run.stderr[-2000:])
        else:
            events=[json.loads(l) for l in run.stdout.splitlines()];end=events[-1]
            report=dict(class_index=item['class_index'],**end,
                        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),log=str(log))
        reports.append(report);print(json.dumps(report),flush=True)
        (BASE/'batch_checkpoint.json').write_text(json.dumps(dict(reports=reports))+'\n')
    result=dict(status='BOUNDED_BATCH_FINISHED',reports=reports,
                complete_Z_exclusions=sum(r.get('Z_excluded',False) for r in reports),
                limit_classes=sum(r['status']=='LIMIT' for r in reports),
                canonical_Z_leaves=sum(r.get('leaves',0) for r in reports),
                seconds=time.monotonic()-start,
                scope='Necessary rank-two pure-Z distance, isotropy and minimum stabilizer conditions; '
                      'complete zero-leaf runs exclude C, leaves need full X extension tests.',
                binary_sha256=hashlib.sha256((BASE/'enumerator').read_bytes()).hexdigest())
    (BASE/'batch.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='reports'}))
if __name__=='__main__':main()
