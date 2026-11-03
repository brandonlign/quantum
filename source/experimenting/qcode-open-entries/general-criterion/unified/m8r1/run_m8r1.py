"""[[17,3,6]] branch m=8,r=1: canonical Z enumeration (C, N=8192, d=6) + exact fixed-Z solver.

Per suffix class: EXCLUDED only if Z enumeration COMPLETE and every leaf is
excluded by z_condition or solve_extension; FOUND requires a verified map;
anything bounded is LIMIT. Minimum stabilizer weight 8 is imposed (branch m=8).
Usage: run_m8r1.py START END [--seconds S] [--max-nodes K] [--solver-nodes K]
"""
import argparse,json,subprocess,sys,tempfile,time,hashlib
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from model import ExtensionModel,nullspace
from signed_extension import solve_extension
from run_benchmark import z_condition
SRC=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
N_,K_,D_,M_=17,3,6,8

BACKENDS = {
    'plain': ('z_enum8192', []),
    'gamma': ('z_enum_xdom', ['--gamma']),
    'gamma_xdom': ('z_enum_xdom', ['--gamma', '--xdom']),
    'rank': ('z_enum_rank', ['--gamma']),
    'rank_xdom': ('z_enum_rank', ['--gamma', '--xdom']),
    'fixed_h5': ('z_enum_fixed_radical', ['--gamma', '--xdom', '--radical-rank', '3']),
    'avoid_h5': ('z_enum_avoid_radical', ['--gamma', '--avoid-radical']),
    'symmetry_h5': ('z_enum_symmetry', ['--gamma', '--avoid-radical', '--first-orbits']),
}

def leaves(distances,seconds,max_nodes,stats,*,Gamma,backend,first_representatives=None):
    binary,flags=BACKENDS[backend]
    if '--first-orbits' in flags and first_representatives is None:
        raise ValueError('first-orbits backend requires verified representatives')
    cmd=[str(HERE/binary),*flags]
    if seconds is not None:cmd+=['--seconds',str(seconds)]
    if max_nodes is not None:cmd+=['--max-nodes',str(max_nodes)]
    with tempfile.TemporaryFile('w+') as err:
        p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1)
        payload=f'{M_} {len(distances)} {D_}\n'+' '.join(map(str,distances))+'\n'
        if '--gamma' in flags:payload+=' '.join(map(str,Gamma))+'\n'
        if '--first-orbits' in flags:
            payload+=str(len(first_representatives))+'\n'+' '.join(map(str,first_representatives))+'\n'
        p.stdin.write(payload);p.stdin.close()
        n=0
        try:
            for line in p.stdout:
                f=line.split();assert f[0]=='LEAF' and len(f)==M_;n+=1
                yield tuple(map(int,f[1:]))
            p.wait();err.seek(0);s=json.loads(err.read().strip());stats.update(s)
            if p.returncode or s.get('status')=='ERROR':
                raise RuntimeError(f'enumerator failed: {s}')
            assert s['status'] in('COMPLETE','LIMIT') and s['leaves']==n
        finally:
            if p.poll() is None:p.terminate();p.wait()

def run_class(C,seconds=None,max_nodes=None,solver_nodes=None,backend='plain'):
    t0=time.time()
    model=ExtensionModel(N_,K_,D_,M_,[((1<<M_)-1)<<M_],C)
    assert model.C_distance>=6 and model.b==13
    stats={};tally=Counter();found=None;solver_limits=0
    for xs in leaves(model.distances,seconds,max_nodes,stats,Gamma=model.Gamma,backend=backend):
        if not z_condition(xs,model.Gamma,model.b):tally['Z']+=1;continue
        r=solve_extension(model,list(xs),node_limit=solver_nodes,minimum_stabilizer_weight=8)
        tally[r['status']+':'+r.get('reason','')]+=1
        if r['status']=='LIMIT':solver_limits+=1
        if r['status']=='FOUND':found=dict(xs=xs,**{k:r[k] for k in['Y','model_columns']});break
    if found:status='FOUND'
    elif stats.get('status')=='COMPLETE' and solver_limits==0:status='EXCLUDED'
    else:status='LIMIT'
    return dict(status=status,backend=backend,z=stats,tally=dict(tally),found=found,hull=model.hull,seconds=round(time.time()-t0,2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('start',type=int);ap.add_argument('end',type=int)
    ap.add_argument('--seconds',type=float);ap.add_argument('--max-nodes',type=int);ap.add_argument('--solver-nodes',type=int)
    ap.add_argument('--only',help='comma list of class indices')
    ap.add_argument('--backend',choices=['plain','gamma','gamma_xdom','rank','rank_xdom'],default='plain')
    ap.add_argument('--output',type=Path,help='JSON checkpoint after each completed class')
    a=ap.parse_args();S=json.loads(SRC.read_text())['classes']
    ids=[int(x) for x in a.only.split(',')] if a.only else range(a.start,min(a.end,len(S)))
    records=[]
    binary,flags=BACKENDS[a.backend]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [SRC,HERE/'run_m8r1.py',HERE/binary]}
    for i in ids:
        r=run_class(S[i]['representative']['generators'],a.seconds,a.max_nodes,a.solver_nodes,a.backend)
        record=dict(cls=i,**r);records.append(record)
        print(json.dumps(record),flush=True)
        if a.output:
            a.output.write_text(json.dumps(dict(records=records,sha256=hashes,
                scope='EXCLUDED requires complete Z enumeration and no unresolved fixed-Z solver'),indent=2)+'\n')
