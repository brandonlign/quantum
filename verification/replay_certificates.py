#!/usr/bin/env python3
"""Replay the compact certificate bundle using only Python's standard library.

This checks the enumerator reductions, not classification or extension searches.
Run without -O, from any working directory. No original project files are needed.
"""
from math import comb
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def kraw(n, j, w):
    return sum((-1)**h * 3**(j-h) * comb(w, h) * comb(n-w, j-h)
               for h in range(max(0, j-(n-w)), min(j, w)+1))


def parity16():
    n, size, modulus = 16, 2048, 4096
    low = [-740, -265, -81, -19, -3]
    shadow = [80, 10, 0, -2]
    coefficients = []
    for w in range(n+1):
        c = low[0]
        c += sum(low[j]*(kraw(n, j, w)-size*int(w == j))
                 for j in range(1, 5))
        c += sum(v*(-1)**w*kraw(n, j, w) for j, v in enumerate(shadow))
        coefficients.append(c % modulus)
    assert coefficients == [size*int(1 <= w <= 4) for w in range(n+1)]
    assert all((-size*v) % modulus == 0 for v in shadow)
    assert (low[0]*size) % modulus == 0
    print('PASS sixteen-qubit parity: A1+A2+A3+A4 is even', flush=True)


def profile17():
    n, size = 17, 16384
    residual, rhs = [0]*18, 0
    def add(mult, row, value=0):
        nonlocal rhs
        residual[:] = [a+mult*b for a, b in zip(residual, row)]
        rhs += mult*value
    add(-1835008, [int(w == 0) for w in range(18)], 1)
    add(116, [1]*18, size)
    for j, mult in enumerate([60, 25, 11, 3, 1], 1):
        add(mult, [kraw(n, j, w)-size*int(j == w) for w in range(18)])
    for j, mult in enumerate([-851968, -278528, -507904, -147456, -180224], 1):
        add(mult, [int(w == j) for w in range(18)])
    for j, mult in enumerate([88, 32, 18, 4, 2]):
        assert mult >= 0
        add(mult, [-(-1)**w*kraw(n, j, w) for w in range(18)])
    assert residual == [32768*int(w in (6, 7)) for w in range(18)]
    assert rhs == 65536
    print('PASS seventeen-qubit profile: A6+A7 <= 2 under purity', flush=True)


def main():
    assert __debug__, 'Assertions must remain enabled; run without -O'
    total = 0
    for group, checker, expected in [
        ('14', 'check_certificates14.py', 11),
        ('16', 'check_certificates.py', 6),
        ('17', 'check_certificates.py', 7),
        ('17-frame-counts', 'check_certificates.py', 13),
    ]:
        paths = sorted((HERE/'certificates'/group).glob('*.json'))
        assert len(paths) == expected, (group, len(paths), expected)
        subprocess.run([sys.executable, '-B', str(HERE/checker), *map(str, paths)], check=True)
        total += len(paths)
    for script in ['low_weight14.py', 'parity14.py']:
        subprocess.run([sys.executable, '-B', str(HERE/script)], check=True)
    parity16()
    profile17()
    print(f'PASS: {total} exact certificates and four ordinary identities', flush=True)


if __name__ == '__main__':
    main()
