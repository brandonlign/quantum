import json, time, sys
import numpy as np
from sympl import setup
from zstage import leaves
cells = json.load(open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json'))['7,4']
for i in (0, 1):
    dmin, syn, Gam, b = setup(cells[i], 7)
    t = time.time(); c = 0
    for xs in leaves(dmin, 9):
        c += 1
        if time.time() - t > 60: break
    print(i, 'leaves in 60s', c, 'rank Gamma', np.linalg.matrix_rank(Gam.astype(float)), 'b', b, 'deep cosets', int((dmin >= 4).sum()), flush=True)
