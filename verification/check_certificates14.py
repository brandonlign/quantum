#!/usr/bin/env python3
"""Standalone checker for the Section 5 Farkas certificates (Appendix A).

Each file in verification/certificates/ specifies one linear system in the
format of Appendix A of the paper and an integer Farkas certificate for it.
This script rebuilds every system from scratch:

  * coset weight patterns are recomputed from the subgroup generators;
  * every constraint row is rebuilt from its label (Appendix A.2);
  * coset-relation rows are checked against the recomputed patterns;
  * the residual  c = sum(lambda*G) + sum(nu*E)  is computed in every
    column and must be >= 0, while  sum(nu*rhs)  must be < 0.

Only the Python standard library and exact integer arithmetic are used.
It shares no code with the other verification scripts.

Usage: python3 verification/check_certificates.py [FILE.json ...]
"""
import json
import sys
from itertools import product
from math import comb, prod
from pathlib import Path

N_STAB = 2048  # |S| = 2^11


def kraw(n, a, i):
    """K_a^n(i) = [u^a] (1+3u)^(n-i) (1-u)^i."""
    return sum((-1) ** h * 3 ** (a - h) * comb(i, h) * comb(n - i, a - h)
               for h in range(max(0, a - (n - i)), min(a, i) + 1))


def parse_pauli(word):
    x = z = 0
    for k, ch in enumerate(word):
        if ch in "XY":
            x |= 1 << k
        if ch in "ZY":
            z |= 1 << k
        assert ch in "IXYZ", word
    return x, z


def coset_patterns(m, generators):
    """Weight histograms of the R-cosets in the normalizer of R in V_m."""
    gens = [parse_pauli(g) for g in generators]
    assert all(len(g) == m for g in generators)
    symp = lambda u, v: ((u[0] & v[1]).bit_count() + (u[1] & v[0]).bit_count()) & 1
    assert all(symp(g, h) == 0 for g in gens for h in gens), "R not isotropic"
    R = {(0, 0)}
    for g in gens:
        R |= {(x ^ g[0], z ^ g[1]) for x, z in R}
    normal = [(x, z) for x in range(1 << m) for z in range(1 << m)
              if all(symp((x, z), g) == 0 for g in gens)]
    seen, patterns = set(), {}
    for v in normal:
        if v in seen:
            continue
        coset = {(v[0] ^ r[0], v[1] ^ r[1]) for r in R}
        seen |= coset
        hist = [0] * (m + 1)
        for x, z in coset:
            hist[(x | z).bit_count()] += 1
        patterns[tuple(hist)] = patterns.get(tuple(hist), 0) + 1
    assert sum(patterns.values()) * len(R) == len(normal)
    return patterns


def check(path):
    spec = json.loads(Path(path).read_text())
    blocks, f = spec["blocks"], spec["suffix_sites"]
    assert sum(b["sites"] for b in blocks) + f == 14
    sizes = [b["sites"] for b in blocks] + [f]

    # Variables: one coset-pattern choice per block, plus a suffix weight.
    block_patterns = []
    for b in blocks:
        m = b["sites"]
        pats = coset_patterns(m, b["generators"]) if b["generators"] else None
        if b["variables"] == "coset_patterns":
            assert pats is not None
            want = {tuple(p): c for p, c in b["patterns"]}
            assert pats == want, f"coset patterns differ in {path}"
            block_patterns.append(sorted(pats))
        else:  # "cells": the variables are the split counts H themselves
            assert b["variables"] == "cells"
            block_patterns.append([tuple(int(a == k) for a in range(m + 1))
                                   for k in range(m + 1)])
        b["_patterns"] = pats
    columns = list(product(*[range(len(p)) for p in block_patterns], range(f + 1)))
    cells = list(product(*[range(n + 1) for n in sizes]))
    assert len(columns) == spec["expected"]["columns"]

    # Per-column factors: H_c(y) = prod P_i[c_i] * [c_f == b];
    # T_c and SH_c are the tensor-Krawtchouk transforms of H.
    def block_vecs(i, P, signed):
        n = sizes[i]
        s = [(-1) ** a if signed else 1 for a in range(n + 1)]
        return [sum(kraw(n, a, j) * s[j] * P[j] for j in range(n + 1))
                for a in range(n + 1)]
    tvec = [[block_vecs(i, P, False) for P in pl] for i, pl in enumerate(block_patterns)]
    svec = [[block_vecs(i, P, True) for P in pl] for i, pl in enumerate(block_patterns)]

    def H(c, col):
        return prod(block_patterns[i][col[i]][c[i]] for i in range(len(blocks))) * (c[-1] == col[-1])

    def T(c, col, signed=False):
        v = svec if signed else tvec
        sfx = kraw(f, c[-1], col[-1]) * ((-1) ** col[-1] if signed else 1)
        return prod(v[i][col[i]][c[i]] for i in range(len(blocks))) * sfx

    def cell(row):
        c = tuple(row["cell"])
        assert len(c) == len(sizes) and all(0 <= c[i] <= sizes[i] for i in range(len(c)))
        return c

    # Build each labelled row as (coefficients over columns, rhs).
    def build(row):
        t = row["type"]
        if t == "unit":                       # H_0 = 1
            return [H((0,) * len(sizes), col) for col in columns], 1
        if t == "total":                      # sum_c H_c = 2048
            return [sum(H(c, col) for c in cells) for col in columns], N_STAB
        if t == "low":                        # T_c - 2048 H_c = 0, 1 <= |c| <= 4
            c = cell(row)
            assert 1 <= sum(c) <= 4
            return [T(c, col) - N_STAB * H(c, col) for col in columns], 0
        if t == "zero":                       # H_c = 0 for |c| in {1,3} (Lemma 1)
            c = cell(row)
            assert sum(c) in (1, 3)
            return [H(c, col) for col in columns], 0
        if t == "weight":                     # sum_{|c|=w} H_c = A_w
            w, value = row["w"], row["value"]
            assert (w, value) in {(1, 0), (2, 0), (3, 0), (4, 3)}
            return [sum(H(c, col) for c in cells if sum(c) == w) for col in columns], value
        if t == "incl":                       # 2048 H_c - T_c <= 0
            c = cell(row)
            return [N_STAB * H(c, col) - T(c, col) for col in columns], 0
        if t == "cnonneg":                    # -T_c <= 0
            c = cell(row)
            return [-T(c, col) for col in columns], 0
        if t == "shadow":                     # -SH_c <= 0
            c = cell(row)
            return [-T(c, col, signed=True) for col in columns], 0
        if t in ("coset_eq", "coset_le"):     # sum_a w_a H_{a,b} (= or <=) 0 on block 0
            assert len(blocks) == 1
            weights, b = row["weights"], row["suffix_weight"]
            pats = blocks[0]["_patterns"]
            for P in pats:                    # implied by every coset pattern
                val = sum(w * p for w, p in zip(weights, P))
                assert (val == 0) if t == "coset_eq" else (val <= 0), (t, P)
            return [sum(w * H((a, b), col) for a, w in enumerate(weights))
                    for col in columns], 0
        raise ValueError(t)

    eq_types = {"unit", "total", "low", "zero", "weight", "coset_eq"}
    le_types = {"incl", "cnonneg", "shadow", "coset_le"}
    residual = [0] * len(columns)
    nu_rhs = 0
    for row in spec["certificate"]:
        mult = row["multiplier"]
        assert isinstance(mult, int) and mult != 0
        if row["type"] in le_types:
            assert mult > 0, "inequality multipliers must be nonnegative"
        else:
            assert row["type"] in eq_types
        coeffs, rhs = build(row)
        for k, v in enumerate(coeffs):
            residual[k] += mult * v
        nu_rhs += mult * rhs if row["type"] in eq_types else 0

    assert all(r >= 0 for r in residual), f"negative residual in {path}"
    assert nu_rhs < 0
    exp = spec["expected"]
    positive = sum(r > 0 for r in residual)
    assert nu_rhs == exp["nu_rhs"] and positive == exp["positive_residuals"]
    print(f"PASS {Path(path).name}: {len(columns)} columns, "
          f"{len(spec['certificate'])} multiplier rows, "
          f"{positive} positive residuals, nu.rhs = {nu_rhs}")


if __name__ == "__main__":
    files = sys.argv[1:] or sorted(
        (Path(__file__).resolve().parent / "certificates").glob("*.json"))
    assert files, "no certificate files found"
    for p in files:
        check(p)
    print(f"ALL {len(files)} SECTION 5 CERTIFICATES PASS")
