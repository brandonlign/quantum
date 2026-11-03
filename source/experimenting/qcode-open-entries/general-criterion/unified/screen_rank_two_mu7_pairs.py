"""Bounded rank-two screen with minimum stabilizer7; preserves original mu8 files."""
import argparse
import hashlib
import json
from pathlib import Path
from model import lift_coordinate
import subprocess

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=10)
    ap.add_argument('--classes');a=ap.parse_args()
    binary=HERE/'rank2_mu7/enumerator_pair_domains';out=HERE/'rank2_mu7_pairs';out.mkdir(exist_ok=True)
    assert json.loads((HERE/'rank2_mu7/pair_domains_verification.json').read_text())['status']=='PASS'
    items=json.loads((HERE/'rank2_z_inputs.json').read_text())['classes']
    selected=None if a.classes is None else set(map(int,a.classes.split(',')))
    orbit_path=HERE/'rank2_z/aut_orbits.json'
    orbit_info=json.loads(orbit_path.read_text());assert orbit_info['status']=='PASS'
    orbit_info={r['class_index']:r for r in orbit_info['classes']}
    records=[]
    for item in items:
        index=item['class_index']
        if selected is not None and index not in selected:continue
        source=Path(item['table']);data=json.loads(source.read_text())
        assert data['m']==8 and data['p']==6
        rows=[f'4096 6 100000000 {int(a.seconds*1000)} 10',' '.join(map(str,data['needs']))]
        for label,(form,delta,lam) in enumerate(zip(data['form_rows'],data['coset_distances'],data['suffix_weights'])):
            flags=(int(delta>=4 and lam>=5)|2*int(delta>=2 and lam>=3)) if label else 0
            mask=sum(1<<c for c in range(9) if delta>=6-c and lam>=7-c)
            if label:assert bool(flags&1)==bool(mask&4) and bool(flags&2)==bool(mask&16)
            rows.append(f'{flags} {form} {mask}')
        auto=orbit_info[index]
        assert auto['table_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
        high={i for i,(delta,lam) in enumerate(zip(data['coset_distances'],data['suffix_weights'])) if i and delta>=4 and lam>=5}
        remaining=set(high);first=[]
        while remaining:
            seed=min(remaining);orbit={seed};front=[seed]
            while front:
                v=front.pop()
                for T in auto['verified_generator_coefficient_maps']:
                    w=lift_coordinate(v,T);assert w in high
                    assert data['coset_distances'][v]==data['coset_distances'][w] and data['suffix_weights'][v]==data['suffix_weights'][w]
                    if w not in orbit:orbit.add(w);front.append(w)
            assert orbit<=remaining;remaining-=orbit;first.append(seed)
        rows.extend([str(len(first)),' '.join(map(str,first))])
        path=out/f'class_{index}.txt';path.write_text('\n'.join(rows)+'\n')
        run=subprocess.run([str(binary),str(path)],text=True,capture_output=True,check=True)
        log=out/f'class_{index}.jsonl';log.write_text(run.stdout)
        stats=json.loads(run.stdout.splitlines()[-1]);assert stats['status'] in ['COMPLETE','LIMIT']
        records.append(dict(class_index=index,first_orbits=len(first),**stats,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),log=str(log)))
        result=dict(status='BOUNDED_MU7_PAIR_SCREEN',minimum_stabilizer=7,records=records,
                    complete_exclusions=sum(r['status']=='COMPLETE' and r['leaves']==0 for r in records),
                    sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            [Path(__file__),orbit_path,binary,HERE/'rank2_mu7/enumerator_pair_domains.cpp']},
                    scope='Singleton and pair-X sector necessary domains; LIMIT/surviving leaves do not exclude; original minimum8 result unchanged')
        (out/'screen.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(records[-1]),flush=True)


if __name__=='__main__':main()
