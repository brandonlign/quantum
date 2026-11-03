import sys, json
sys.path.insert(0,'split'); sys.path.insert(0,'general-criterion/unified/frame_selection')
from certs import with_row, pure, n,k,d   # reuse helper (re-runs certs on import is avoided below)
