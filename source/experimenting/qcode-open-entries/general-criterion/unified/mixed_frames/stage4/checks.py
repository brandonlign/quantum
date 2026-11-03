"""Serial, bounded subprocess calls; every input and output is Stage-4-owned."""
import json
import subprocess
from .common import HERE,BINARY,input_text


def cpp(model,frame,*,mu=0,counts=None,w6=6,w7=7,prefix=(),depth=None,
        nodes=None,leaves=None,seconds=None,rank=True,domains=True,binary=BINARY,
        canonical=True,symmetry=None):
    path=HERE/'validation_input.txt'
    path.write_text(input_text(model,frame,mu,rank,domains,counts,w6,w7))
    cmd=[str(binary),str(path)]
    options=[('prefix',','.join(map(str,prefix)) if prefix else None),('depth',depth),
        ('nodes',nodes),('leaves',leaves),('seconds',seconds)]
    if binary==BINARY:options += [('canonical',int(canonical)),('symmetry',symmetry)]
    for k,v in options:
        if v is not None:cmd+=['--'+k,str(v)]
    proc=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=15)
    rows=[json.loads(line) for line in proc.stdout.splitlines()]
    assert rows and rows[-1]['type']=='summary' and sum(r['type']=='summary' for r in rows)==1
    assert all(r['type']=='leaf' for r in rows[:-1])
    return dict(**rows[-1],leaves=[r['X'] for r in rows[:-1]])
