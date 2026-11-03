import json, sys, subprocess, os
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
from full import syndromes
from phi_input import make_input
# census layout is interleaved symbols; cells.json layout is X bits | Z bits<<n
def deinter(v, n):
    x = z = 0
    for i in range(n):
        x |= ((v >> (2*i)) & 1) << i; z |= ((v >> (2*i+1)) & 1) << i
    return x | (z << n)
cells = json.load(open('cells11.json'))['11,12']
for i, g in enumerate(cells):
    gens = [deinter(v, 11) for v in g]
    dmin, syn, r = syndromes(gens, 11)
    print('class', i, 'W/C dim', r, flush=True)
    assert r == 10
    out = subprocess.run(['./phisearch2', '0'], input=make_input(dmin, 6), capture_output=True, text=True).stdout.strip()
    print(6, i, out, flush=True)
