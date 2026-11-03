import sys
import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']
exec(open('control3.py').read().split("ok, st = search")[0])
from sat_check import own_dmin
from sat_sym import build
from pysat.solvers import Cadical153
d2, r2 = own_dmin(Cg, f)
with Cadical153(bootstrap_with=build(d2, r2, m)) as S:
    print('symmetry-broken SAT control on [[16,4,5]] (r=%d, m=%d):' % (r2, m), S.solve(), '(must be True)')
