import json, sys, time, collections
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general')
from full import syndromes
from hierarchy import screen
cells = json.load(open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json'))
n = int(sys.argv[1]); m = int(sys.argv[2]); key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}[m - (n - 16)]
c = collections.Counter(); surv = []; t = time.time()
for i, g in enumerate(cells[key]):
    try: res, nodes = screen(syndromes(g, n - m)[0], m, node_limit=2 * 10**6)
    except TimeoutError: res = 'LIMIT'
    c[res] += 1
    if res in ('PASS', 'LIMIT'): surv.append((i, res))
print(f'n={n} m={m} cell {key}: {dict(c)} survivors {surv} {time.time()-t:.0f}s', flush=True)
