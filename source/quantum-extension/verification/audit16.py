"""Independent evaluation audit; imports no project proof modules.

Checks all saved suffix codebooks and terminal-log coverage. Reconstructs coset
distances by BFS on the quotient Cayley graph, then builds physical prefix
constraints directly for a bounded fresh replay of the compiled searches.
This does not reclassify the suffix cells or rerun every long search.
"""
import collections
import hashlib
import json
import pathlib
import re
import subprocess
import tempfile
import time

PACKAGE = pathlib.Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
OUT = PACKAGE / 'audit'
PKG = PACKAGE / 'verification/release16'
CELLS = json.loads((PKG / 'census/cells.json').read_text())

def search_binary():
    """Build the retained production C search; only its inputs are independent."""
    source = pathlib.Path(__file__).resolve().parent / 'native/phisearch2.c'
    folder = tempfile.TemporaryDirectory(prefix='quantum-extension-16-')
    binary = pathlib.Path(folder.name) / 'phisearch2'
    subprocess.run(['cc', '-O3', '-o', str(binary), str(source)], check=True)
    return folder, binary

def nullspace(rows, width):
    rows = list(rows)
    pivots = []
    rank = 0
    for col in range(width):
        pivot = next((j for j in range(rank, len(rows)) if rows[j] >> col & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for j in range(len(rows)):
            if j != rank and rows[j] >> col & 1:
                rows[j] ^= rows[rank]
        pivots.append(col)
        rank += 1
    result = []
    for col in range(width):
        if col in pivots:
            continue
        v = 1 << col
        for j, p in enumerate(pivots):
            if rows[j] >> col & 1:
                v |= 1 << p
        result.append(v)
    assert all((h & g).bit_count() % 2 == 0 for h in result for g in rows)
    return rank, result

def quotient_distances(gens, f):
    rank, checks = nullspace(gens, 2*f)
    assert rank == len(gens)
    def syndrome(w):
        return sum(((w & h).bit_count() % 2) << i for i, h in enumerate(checks))
    edges = set()
    for i in range(f):
        a, b = syndrome(1 << i), syndrome(1 << (i+f))
        edges.update((a, b, a ^ b))
    dist = [-1] * (1 << len(checks))
    dist[0] = 0
    queue = collections.deque([0])
    while queue:
        s = queue.popleft()
        for edge in edges:
            t = s ^ edge
            if dist[t] < 0:
                dist[t] = dist[s] + 1
                queue.append(t)
    assert min(dist) == 0 and dist.count(0) == 1
    return dist

def physical_constraints(m):
    mask = (1 << m) - 1
    short = (1 << (m-1)) - 1
    constraints = []
    for c in range(1, 1 << (2*m-2)):
        z = c & short
        x = c >> (m-1)
        x |= (x.bit_count() % 2) << (m-1)
        cost = min((x | z).bit_count(), (x | (z ^ mask)).bit_count())
        if cost < 5:
            constraints.append((c, 5-cost))
    return constraints

def log_records(name, pattern):
    text = (PKG / 'search/logs' / name).read_text()
    return [(int(m), int(i), status) for m, i, status in re.findall(pattern, text)]

def main():
    assert __debug__, 'Run without -O'
    OUT.mkdir(exist_ok=True)
    build, binary = search_binary()
    start = time.monotonic()
    report = {'scope': __doc__, 'cells': {}, 'fresh_searches': []}
    for m, key, expected in [(6,'10,10',37),(7,'9,8',2298),(8,'8,6',1579),(9,'7,4',43)]:
        codes = CELLS[key]
        assert len(codes) == expected
        f = 16-m
        rank_expected = 22-2*m
        for gens in codes:
            assert len(gens) == rank_expected
            assert nullspace(gens,2*f)[0] == rank_expected
            span = [0]
            for g in gens:
                span += [w ^ g for w in span]
            assert len(set(span)) == 1 << rank_expected
            assert min(((w & ((1<<f)-1)) | (w >> f)).bit_count() for w in span[1:]) >= 5
        report['cells'][str(m)] = {'classes':expected,'all_dimensions_and_distances_pass':True}

    dfs = log_records('full2_run.log',r'm=(\d+) class (\d+): (excluded|SURVIVES)')
    dfs += log_records('run_c78.log',r'(?m)^(\d+) (\d+) (NONE|LIMIT|FOUND)')
    assert (PKG/'search/logs/m8_425_rerun.log').read_text().startswith('NONE nodes=23269942')
    dfs = [(m,i,'NONE' if m==8 and i==425 and s=='LIMIT' else s) for m,i,s in dfs]
    dfs += log_records('run_c9.log',r'(?m)^(\d+) (\d+) (NONE|LIMIT|FOUND)')
    second = log_records('sat_run.log',r'm=(\d+) class (\d+): (UNSAT|SAT)')
    second += log_records('sat_sym_78.log',r'm=(\d+) class (\d+): (UNSAT|SAT)')
    for name in ['run_m9_alt_0-1.log','run_m9_alt_2-42.log']:
        second += log_records(name,r'(?m)^(\d+) (\d+) (NONE|LIMIT|FOUND)')
    for label, records in [('primary',dfs),('second',second)]:
        for m in range(6,10):
            cell = [(i,s) for mm,i,s in records if mm==m]
            expected = report['cells'][str(m)]['classes']
            assert all(s in ['NONE','excluded','UNSAT'] for i,s in cell), (label,m)
            ids = sorted({i for i,s in cell})
            assert ids==list(range(expected)), (label,m,len(ids),expected)
            report['cells'][str(m)][label+'_logged_records']=len(cell)
        report[label+'_terminal_log_coverage'] = 3957
    print('All 3957 codebook dimensions/distances and both 3957-case log coverages PASS',flush=True)

    constraints = {m:physical_constraints(m) for m in range(6,10)}
    selected = [(6,i) for i in range(37)] + [(7,i) for i in [0,979,2297]] + [(8,i) for i in [0,1578]]
    for m,i in selected:
        key = {6:'10,10',7:'9,8',8:'8,6',9:'7,4'}[m]
        distances = quotient_distances(CELLS[key][i],16-m)
        assert len(distances)==1024
        rows = constraints[m]
        data = f'{m} 1024\n'+' '.join(map(str,distances))+f'\n{len(rows)}\n'+''.join(f'{c} {t}\n' for c,t in rows)
        run = subprocess.run([str(binary)],input=data,text=True,capture_output=True,timeout=120)
        assert run.returncode==0 and re.fullmatch(r'NONE nodes=\d+\n?',run.stdout), (m,i,run.stdout,run.stderr)
        report['fresh_searches'].append({'m':m,'class':i,'status':run.stdout.strip(),'input_sha256':hashlib.sha256(data.encode()).hexdigest()})
        print(m,i,run.stdout.strip(),flush=True)
    report['elapsed_seconds'] = time.monotonic()-start
    report['status']='PASS'
    (OUT/'audit16.json').write_text(json.dumps(report,indent=2)+'\n')
    build.cleanup()
    print('INDEPENDENT INPUT AUDIT PASS',flush=True)

if __name__=='__main__':
    main()
