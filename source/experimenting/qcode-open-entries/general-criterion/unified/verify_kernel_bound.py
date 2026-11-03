"""Finite checks for the binary projection and balanced rank-two kernel bounds."""
from collections import Counter, deque
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from kernel_bound import griesmer_max_dimension, kernel_screen
from model import (coordinates, dual, gram, lift_coordinate, nullspace, reduced_basis,
                   recover_extension, span, weight)
from verify_unified import pauli

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def verify():
    binary_count = 0
    for length in range(1, 7):
        seen = {()}
        queue = deque([()])
        while queue:
            B = queue.popleft()
            if B:
                distance = min(v.bit_count() for v in span(B) if v)
                assert len(B) <= griesmer_max_dimension(length, distance)
            binary_count += 1
            for v in range(1, 1 << length):
                next_basis = tuple(reduced_basis(list(B) + [v]))
                if next_basis not in seen:
                    seen.add(next_basis)
                    queue.append(next_basis)
    S = [pauli(s) for s in ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']]
    five, cols = recover_extension(S, 5, 1, 3, 4)
    Kcoords = nullspace([sum(((c >> j) & 1) << i for i, c in enumerate(cols))
                        for j in range(five.b)], five.a)
    Kphys = [lift_coordinate(v, five.E) for v in Kcoords]
    assert len(Kphys) == 2
    # The control's two full-support local words are not initially X^4,Z^4.
    # Normalize both simultaneously by the local symplectic coordinate change.
    p0, q0 = five.R
    def normalize(v):
        x = z = 0
        for i in range(4):
            vx, vz = (v >> i) & 1, (v >> (4 + i)) & 1
            px, pz = (p0 >> i) & 1, (p0 >> (4 + i)) & 1
            qx, qz = (q0 >> i) & 1, (q0 >> (4 + i)) & 1
            assert (px * qz + pz * qx) % 2 == 1
            x |= ((px * vz + pz * vx) % 2) << i
            z |= ((qx * vz + qz * vx) % 2) << i
        return x | (z << 4)
    assert {normalize(v) for v in five.R} == {15, 15 << 4}
    Kphys = list(map(normalize, Kphys))
    rank_two = []
    for s in range(1, 5):
        m = 4 * s

        def repeat(v):
            return sum(((v >> i) & 1) << (i * s + j) for i in range(4) for j in range(s)) | (
                sum(((v >> (4 + i)) & 1) << (i * s + j) for i in range(4) for j in range(s)) << m)

        R = [(1 << m) - 1, ((1 << m) - 1) << m]
        K = list(map(repeat, Kphys))
        L = R + K
        assert len(reduced_basis(L)) == 4
        spectrum = Counter(weight(v, m) for v in span(L))
        assert spectrum == {0: 1, 3 * s: 12, 4 * s: 3}
        t = len(reduced_basis(gram(K, m)))
        assert t == (2 if s % 2 else 0)
        for v in span(K):
            if not v:
                continue
            counts = Counter(((v >> i) & 1) + 2 * ((v >> (m + i)) & 1) for i in range(m))
            assert sorted(counts.values()) == [s] * 4
        N = 16
        formula = Fraction(6 * m * m, N) - Fraction(3 * m, 2)
        k2 = lambda w: sum((-1)**j * 3**(2-j) * comb(w, j) * comb(m-w, 2-j)
                           for j in range(max(0, 2 - (m - w)), min(2, w) + 1))
        direct = Fraction(sum(k2(w) * count for w, count in spectrum.items()), N)
        assert formula == direct and direct >= 0
        rank_two.append(dict(m=m, distance=3*s, kernel_dimension=2,
                             kernel_form_rank=t, dual_weight2_count=str(direct)))
    control = kernel_screen(5, 1, 3, 4, 2, five.hull)
    assert control['status'] == 'UNRESOLVED'
    assert any(v['kernel_dimension'] == 2 and v['kernel_form_rank'] == 2
               for v in control['possibilities'])
    fourteen = kernel_screen(14, 3, 5, 4, 1)
    assert fourteen['allowed_hulls'] == [4]
    target17 = {f'm{m}_r{r}': kernel_screen(17, 3, 6, m, r, minimum_stabilizer_weight=m)
                for m in [6, 7, 8] for r in ([1] if m % 2 else [1, 2])}
    assert target17['m6_r1']['allowed_hulls'] == [3, 5]
    assert target17['m7_r1']['allowed_hulls'] == [1, 3]
    assert target17['m6_r2']['allowed_hulls'] == [4]
    assert target17['m8_r2']['kernel_dimension_bound'] == 3
    assert target17['m8_r2']['allowed_hulls'] == [0, 2, 4]
    assert max(v['kernel_dimension'] for v in target17['m8_r2']['possibilities']) == 2
    nine = kernel_screen(16, 5, 5, 9, 1, minimum_stabilizer_weight=9)
    assert {v['kernel_dimension'] for v in nine['possibilities']} == {6}
    files = [Path(__file__), BASE / 'kernel_bound.py', BASE / 'model.py']
    report = dict(status='PASS', exhaustive_binary_subspaces=binary_count,
                  repeated_balanced_kernel_controls=rank_two,
                  five_qubit_positive_control=control, fourteen_qubit_hull=fourteen,
                  target17_kernel_branches=target17,
                  sixteen_m9_minimum_word=nine,
                  sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  scope='Finite classical-bound checks and explicit balanced local spaces; '
                        'rank-two repeats are local kernel examples, not new global codes. '
                        'Target17 branches remain unclassified and unresolved.')
    (BASE / 'kernel_bound_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('target17_kernel_branches', 'sha256')}))


if __name__ == '__main__':
    verify()
