import sys, os
here = '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search'
sys.path.insert(0, here); sys.path.insert(0, os.getcwd())
cwd = os.getcwd(); os.chdir(here); sys.argv = ['x', 'q1645.txt']; __file__ = here + '/control3.py'
src = open('control3.py').read().split('ok, st = search')[0].replace("import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']", '')
exec(src); os.chdir(cwd)
from zstage import run
print('control [[16,4,5]] m=6 (coisotropic must be > 0):', run(Cg, f, m, first=True))
import numpy as np
from sympl import setup, vec
dm2, syn2, Gam, b = setup(Cg, f)
phi = {}
for x, z in Sp:
    e = prefcoords(x, z); s = syn2(tosuf(x & ~Tm, z & ~Tm))
    assert phi.setdefault(e, s) == s
Zi = [phi[1 << (m + i)] for i in range(m)]
Eta = [phi[(1 << i) | (1 << (m - 1))] for i in range(m - 1)]
G2 = sum(np.outer(vec(Zi[i], b), vec(Eta[i], b)) + np.outer(vec(Eta[i], b), vec(Zi[i], b)) for i in range(m - 1)) % 2
print('identity Gamma = sum(xi eta^T + eta xi^T):', bool((G2 == Gam).all()), ' rank Gamma', int(np.linalg.matrix_rank(Gam.astype(float))))
from sympl import z_coisotropic
print('own xi coisotropic:', z_coisotropic(Zi, Gam, b), ' xi sum zero:', np.bitwise_xor.reduce(Zi) == 0)
