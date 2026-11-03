"""Positive control for phisearch on the existing [[16,4,5]] code (W/C dim 11)."""
import sys, subprocess
import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']
exec(open('control3.py').read().split("ok, st = search")[0])
from phi_input import make_input
import os; PROG = [os.environ.get('PHIPROG', './phisearch2')]
r = subprocess.run(PROG, input=make_input(dmin, m), capture_output=True, text=True)
print('C control on [[16,4,5]]:', r.stdout.strip(), '(must be FOUND)')
