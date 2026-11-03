import json, time, sys
import xstage
from xstage import *
cells = json.load(open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json'))['7,4']
dmin, syn, Gam, b = setup(cells[36], 7)
basis, cost = prefix(9)
for xs in leaves(dmin, 9):
    xi = xs + [int(np.bitwise_xor.reduce(np.array(xs)))]
    t = time.time(); print('U only:', screen_leaf(xs, dmin, Gam, b, 9, cost, basis[:8], basis[8:], pairs=False), f'{time.time()-t:.1f}s', flush=True)
    t = time.time(); print('with pairs:', screen_leaf(xs, dmin, Gam, b, 9, cost, basis[:8], basis[8:]), f'{time.time()-t:.1f}s', flush=True)
