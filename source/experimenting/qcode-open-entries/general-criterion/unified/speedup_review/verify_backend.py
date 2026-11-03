"""Complete small-tree comparisons, physical controls, and limit semantics."""
import itertools
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from model import ExtensionModel, coordinates, lift_coordinate, reduced_basis
from mixed_frames.frames import adapted_frame, mixed_frame
from mixed_frames.stage3.common import input_text, graph_counts, sha
from mixed_frames.stage3.run import complete_leaf
from mixed_frames.fixed_z import validate_solution
from speedup_review.verify_688_symmetry import normal, permute, main as physical_check


def call(text, canonical, binary=None, **options):
    path = HERE / 'validation_input.txt'
    path.write_text(text)
    cmd = [str(binary or HERE/'enumerate_688'), str(path)]
    if binary is None:
        cmd += ['--canonical', str(int(canonical))]
    for k, v in options.items():
        if k == 'prefix': v = ','.join(map(str, v))
        cmd += ['--'+k, str(v)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
    rows = [json.loads(s) for s in result.stdout.splitlines()]
    assert rows[-1]['type'] == 'summary'
    assert sum(r['type'] == 'summary' for r in rows) == 1
    xs = [tuple(r['X']) for r in rows[:-1]]
    assert len(xs) == rows[-1]['leaf_count'] == len(set(xs))
    return rows[-1], set(xs)


def main():
    started = time.monotonic()
    physical_check()
    report = dict(status='PASS', fixtures=[], physical_controls=[], target_checks=[])
    # Uniform suffix data isolate traversal and canonical coverage. The frame
    # is physically R=<Z^6,X^4II>; its S4 action is the same anchor mechanism.
    f = adapted_frame(6, [63<<6, 15])
    assert f.E[:4] == [16<<6, 9<<6, 10<<6, 12<<6]
    for rank in (False, True):
        for domains in (False, True):
            # The disabled symmetry lane must match the frozen backend exactly.
            arrays = [[4,4,3,0,int(rank),int(domains)],f.costs,[6]*16,[0]*4,[20]*16]
            text = '\n'.join(' '.join(map(str,a)) for a in arrays)+'\n'
            a, aa = call(text, False, binary=ROOT/'mixed_frames/stage3/enumerate_z')
            b, bb = call(text, False)
            assert aa == bb and a['nodes'] == b['nodes'] and a['status'] == b['status'] == 'COMPLETE'
            c, cc = call(text, True)
            assert c['status'] == 'COMPLETE'
            normalized = {tuple(normal(xs)) for xs in aa}
            assert normalized == cc and len(aa) == 24*len(cc)
            report['fixtures'].append(dict(rank=rank,domains=domains,plain_leaves=len(aa),canonical_leaves=len(cc),plain_nodes=a['nodes'],canonical_nodes=c['nodes']))
    # Six-site full-plane positive controls exercise complete Y enumeration
    # and the exact-profile final check. S2 anchor sorting is trivial here.
    m = ExtensionModel(6,0,3,4,[15,240],[])
    sf = adapted_frame(4,m.R)
    positive = complete_leaf(m,sf,[3,12],3,None,3,4)
    assert positive['accepted']
    for measured in positive['measurements']:
        if measured['minimum'] < 3: continue
        counts = (measured['n6'],measured['n7'])
        text = input_text(m,sf,3,counts=counts,w6=3,w7=4)
        summary, xs = call(text,True,prefix=[3,12])
        assert summary['status']=='COMPLETE' and xs=={(3,12)}
        done = complete_leaf(m,sf,[3,12],3,counts,3,4)
        assert done['accepted']
    report['physical_controls'].append(dict(n=6,d=3,solutions=len(positive['accepted']),exact_profile=True))
    # Actual eight-site 688 positive graph, including a nontrivial S6 move.
    tf = mixed_frame('688')
    m = ExtensionModel(14,0,2,8,tf.R,[])
    columns = [m.syndrome(1<<(6+i)) for i in range(6)] + [m.syndrome(1<<i) for i in range(6)]
    validate_solution(m,tf,columns)
    wanted = normal(columns[:6])
    transformed = None
    for order in itertools.permutations(range(6)):
        A = [coordinates(permute(w,[*order,6,7]),tf.R+tf.E)>>2 for w in tf.E]
        changed = [lift_coordinate(a,columns) for a in A]
        if changed[:6] == wanted:
            transformed = changed
            break
    assert transformed is not None
    validate_solution(m,tf,transformed)
    measured = graph_counts(m,tf,columns,2,3)
    assert measured == graph_counts(m,tf,transformed,2,3)
    text = input_text(m,tf,2,counts=(measured['n6'],measured['n7']),w6=2,w7=3)
    summary, xs = call(text,True,prefix=wanted)
    assert summary['status']=='COMPLETE' and xs=={tuple(wanted)}
    report['physical_controls'].append(dict(n=14,frame='688',d=2,original=columns[:6],representative=wanted,exact_profile=True))
    # Every actual target hull regime: fixed complete prefixes plus disabled
    # symmetry regression on bounded identical traversal. No exclusion claims.
    manifest = json.loads((ROOT/'mixed_frames/stage3/manifest.json').read_text())
    for index in (0,7,12):
        row = next(r for r in manifest['classes'] if r['key']==f'688_{index}_20')
        text = (ROOT/'mixed_frames/stage3'/row['input']).read_text()
        a,aa = call(text,False,binary=ROOT/'mixed_frames/stage3/enumerate_z',nodes=30)
        b,bb = call(text,False,nodes=30)
        assert aa==bb and a['nodes']==b['nodes'] and a['status']==b['status']
        # Canonical rejection must not report a false exhausted root if capped.
        c,cc = call(text,True,nodes=1)
        assert c['status']=='LIMIT' and not c['complete']
        report['target_checks'].append(dict(index=index,hull=row['hull'],disabled_equal=True,limit=True))
    for kw in (dict(nodes=0),dict(leaves=0),dict(seconds=0)):
        summary,_ = call(input_text(ExtensionModel(6,0,3,4,[15,240],[]),sf),True,**kw)
        assert summary['status']=='LIMIT' and not summary['complete']
    report['seconds'] = time.monotonic()-started
    report['backend_sha256'] = sha(HERE/'enumerate_688')
    report['source_sha256'] = sha(HERE/'enumerate_688.cpp')
    (HERE/'backend_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__ == '__main__': main()
