#!/usr/bin/env bash
# Replay the computational parts of the proof that no [[16,5,5]] qubit stabilizer code exists.
#
#   ./replay.sh            quick replay (about 25 minutes on an Apple M3): certificates, parity
#                          congruence, controls, the full m=6 and m=7 searches, spot checks for m=8, m=9
#   ./replay.sh --census   also re-verify completeness of census/cells.json (needs pynauty;
#                          CENSUS_PYTHON=... selects a different interpreter for this step)
#   ./replay.sh --full     also rerun the full m=8 DFS (~6 core-hours) and the full SAT cross-check
#                          for m=7,8 (~11 core-hours); m=9 is not run here (about 14 core-hours:
#                          cd search && PHIPROG=./phisearch python3 run_c.py 9)
#   Environment: PYTHON (needs numpy, sympy, python-sat), CENSUS_PYTHON, CC.
set -euo pipefail
cd "$(dirname "$0")"
PY=$(command -v "${PYTHON:-python3}" || echo "${PYTHON:-python3}"); case "$PY" in /*) ;; *) PY="$PWD/$PY";; esac
# the census step needs pynauty; set CENSUS_PYTHON if that lives in a different interpreter
PYC=$(command -v "${CENSUS_PYTHON:-$PY}" || echo "${CENSUS_PYTHON:-$PY}"); case "$PYC" in /*) ;; *) PYC="$PWD/$PYC";; esac
CC=${CC:-cc}
CENSUS=0; FULL=0
for a in "$@"; do case "$a" in --census) CENSUS=1;; --full) FULL=1;; *) echo "unknown option $a"; exit 2;; esac; done

echo "== 1. Farkas certificates (Proposition 3; exact integers, standard library only)"
$PY certificates/check_certificates.py

echo "== 2. Parity congruence (Lemma 4)"
$PY parity_congruence.py | tail -1

if [ $CENSUS = 1 ]; then
  echo "== 3. Census: completeness of cells.json and empty cells (Section 5)"
  $PYC census/check_cells.py | tail -6
fi

echo "== 4. Build the C searches"
cd search
for p in phisearch phisearch2 phisearch3; do $CC -O3 -o $p $p.c; done

echo "== 5. Positive controls on the existing [[16,4,5]] code (all must be FOUND/True)"
for p in phisearch phisearch2 phisearch3; do PHIPROG=./$p $PY control_c.py | tail -1; done
$PY sat_control.py | tail -1
$PY sat_sym_control.py | tail -1

echo "== 6. m=6 and m=7: full DFS (no symmetry reduction) and SAT (with Lemma 10)"
PHIPROG=./phisearch2 $PY run_c.py 6 7 | grep survivors
$PY sat_sym.py 6:all | tail -1

echo "== 7. Spot checks m=8 (DFS) and m=9 (DFS with Lemma 10; a few minutes each)"
$PY - <<'EOF'
import json, subprocess
from full import syndromes
from phi_input import make_input
from run_c import cells
for prog, m, key, i in [('./phisearch2', 8, '8,6', 0), ('./phisearch2', 8, '8,6', 1578), ('./phisearch3', 9, '7,4', 42)]:
    dmin, _, r = syndromes(cells[key][i], 16 - m); assert r == 10
    out = subprocess.run([prog], input=make_input(dmin, m), capture_output=True, text=True).stdout.strip()
    print(f'm={m} class {i}: {out}'); assert out.startswith('NONE')
EOF

if [ $FULL = 1 ]; then
  echo "== 8. Full m=8 DFS and full SAT cross-check for m=7,8"
  PHIPROG=./phisearch2 $PY run_c.py 8 | grep survivors
  $PY sat_sym.py 7:all 8:all | tail -1
fi
echo "REPLAY OK"
