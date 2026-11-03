import json, sys, time, collections
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general')
import os
from full import syndromes
from hierarchy import screen
here = '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search'
# positive control: [[16,4,5]] code
os.chdir(here); sys.argv = ['x', 'q1645.txt']
__file__ = here + '/control3.py'
src = open('control3.py').read().split('ok, st = search')[0].replace("import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']", '')
exec(src)
print('control [[16,4,5]] m=6:', screen(dmin, 6), '(must be PASS)', flush=True)
cells = json.load(open('../census/cells.json'))
for (n, key) in [(16, None), (17, None)]:
    for m, key in {16: [(6, '10,10'), (7, '9,8'), (8, '8,6'), (9, '7,4')], 17: [(7, '10,10'), (8, '9,8'), (9, '8,6'), (10, '7,4')]}[n]:
        t = time.time(); c = collections.Counter()
        for g in cells[key]:
            dm = syndromes(g, n - m)[0]
            c[screen(dm, m, node_limit=10**6)[0]] += 1
        print(f'[[{n},{n-11 if n==16 else 6},5]] m={m} cell {key}: {dict(c)}  {time.time()-t:.0f}s', flush=True)
