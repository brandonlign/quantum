"""Bounded rank-two screen with minimum stabilizer6 (weight-eight recentered full-plane frame); preserves original mu8 files."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=2)
    ap.add_argument('--classes');a=ap.parse_args()
    binary=HERE/'rank2_z/enumerator_all_unary';out=HERE/'rank2_mu6';out.mkdir(exist_ok=True)
    assert json.loads((HERE/'rank2_z/all_forward_verification.json').read_text())['status']=='PASS'
    items=json.loads((HERE/'rank2_z_inputs.json').read_text())['classes']
    selected=None if a.classes is None else set(map(int,a.classes.split(',')))
    records=[]
    for item in items:
        index=item['class_index']
        if selected is not None and index not in selected:continue
        source=Path(item['table']);data=json.loads(source.read_text())
        assert data['m']==8 and data['p']==6
        rows=[f'4096 6 100000000 {int(a.seconds*1000)} 10',' '.join(map(str,data['needs']))]
        for label,(form,delta,lam) in enumerate(zip(data['form_rows'],data['coset_distances'],data['suffix_weights'])):
            flags=(int(delta>=4 and lam>=4)|2*int(delta>=2 and lam>=2)) if label else 0
            mask=sum(1<<c for c in range(9) if delta>=6-c and lam>=6-c)
            if label:assert bool(flags&1)==bool(mask&4) and bool(flags&2)==bool(mask&16)
            rows.append(f'{flags} {form} {mask}')
        path=out/f'class_{index}.txt';path.write_text('\n'.join(rows)+'\n')
        run=subprocess.run([str(binary),str(path)],text=True,capture_output=True,check=True)
        log=out/f'class_{index}.jsonl';log.write_text(run.stdout)
        stats=json.loads(run.stdout.splitlines()[-1]);assert stats['status'] in ['COMPLETE','LIMIT']
        records.append(dict(class_index=index,**stats,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),log=str(log)))
        result=dict(status='BOUNDED_MU6_SCREEN',minimum_stabilizer=6,records=records,
                    complete_exclusions=sum(r['status']=='COMPLETE' and r['leaves']==0 for r in records),
                    sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            [Path(__file__),binary,HERE/'rank2_z/enumerator_all_unary.cpp']},
                    scope='Wider local rank-two subproblem; LIMIT/surviving unary leaves do not exclude; original minimum8 result unchanged')
        (out/'screen.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(records[-1]),flush=True)


if __name__=='__main__':main()
