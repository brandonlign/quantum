import sys
import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']
exec(open('control3.py').read().split("ok, st = search")[0])
from sat_check import own_dmin, sat
d2, r2 = own_dmin(Cg, f)
from collections import Counter
assert Counter(d2) == Counter(dmin)
print('SAT control on [[16,4,5]]:', sat(d2, r2, m)[0], '(must be True)')
