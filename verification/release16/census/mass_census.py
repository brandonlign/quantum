#!/usr/bin/env python3
"""Section 7.3: orbit-stabilizer mass check for the length-ten census.

This script recomputes the lengthening census of Section 7 with a different
graph encoding from s7_census.py, and checks in every cell that

    sum_P (m! 6^m / |Aut P|) ext_r(P)  =  sum_D ((m+1)! 6^(m+1) / |Aut D|),

where P runs over parent classes, ext_r(P) counts labelled rank-r
extensions keeping minimum distance 5, and D runs over the classes found.
Each labelled code D is counted exactly once on the left, by its shortening
at the last coordinate. Equality in a cell therefore shows that no class
was lost or counted twice there. It does not use any published count.

It then checks that the 37 classes found are the 37 committed
representatives in s7_37_classes.json, as sets of canonical certificates.

Encoding: a code is a list of generators in F2^(2n); symbol i is bits
(2i, 2i+1). Equivalence is coordinate permutation together with any
permutation of the three nonzero symbols at each coordinate (GL(2,2)).
Graph: codeword vertices (colour 0) and (coordinate, nonzero symbol)
vertices (colour 1); a codeword is joined to the symbols it uses, and the
three symbol vertices of a coordinate form a triangle.

Requires pynauty==2.8.8.1. Runtime: about two minutes.
Usage: python3 verification/s7_mass_check.py [--generators FILE]
"""
import argparse
import json
from fractions import Fraction
from itertools import combinations, product
from math import factorial
from pathlib import Path

from pynauty import Graph, autgrp, certificate

D = 5  # target minimum distance


def span(g):
    out = [0]
    for v in g:
        out += [x ^ v for x in out]
    return out


def wt(v, n):
    return sum(1 for i in range(n) if (v >> (2 * i)) & 3)


def rref(vs):
    piv = {}
    for v in vs:
        x = v
        for b in sorted(piv, reverse=True):
            if (x >> b) & 1:
                x ^= piv[b]
        if x:
            b = x.bit_length() - 1
            for bb in list(piv):
                if (piv[bb] >> b) & 1:
                    piv[bb] ^= x
            piv[b] = x
    return piv


def quotient_labels(gens, n):
    """Label map F2^(2n) -> F2^(2n-k) with kernel <gens>: reduce by RREF,
    then compress the non-pivot bits. Returns function label(v)."""
    piv = rref(gens)
    free = [b for b in range(2 * n) if b not in piv]
    pivs = sorted(piv, reverse=True)

    def label(v):
        x = v
        for b in pivs:
            if (x >> b) & 1:
                x ^= piv[b]
        out = 0
        for i, b in enumerate(free):
            if (x >> b) & 1:
                out |= 1 << i
        return out

    def lift(lab):
        x = 0
        for i, b in enumerate(free):
            if (lab >> i) & 1:
                x |= 1 << b
        return x
    return label, lift, len(free)


def light_errors(n, maxw):
    for w in range(maxw + 1):
        for sites in combinations(range(n), w):
            for syms in product((1, 2, 3), repeat=w):
                yield sum(s << (2 * i) for i, s in zip(sites, syms))


def graph(gens, n):
    words = span(gens)
    W = len(words)
    adj = {v: [] for v in range(W + 3 * n)}
    for r, w in enumerate(words):
        for i in range(n):
            s = (w >> (2 * i)) & 3
            if s:
                adj[r].append(W + 3 * i + s - 1)
    for i in range(n):
        a, b, c = (W + 3 * i + t for t in range(3))
        adj[a] += [b, c]
        adj[b] += [c]
    return Graph(W + 3 * n, directed=False, adjacency_dict=adj,
                 vertex_coloring=[set(range(W)), set(range(W, W + 3 * n))])


def aut_order(g):
    _, s1, s2, _, _ = autgrp(g)
    return round(s1 * 10 ** s2)


def extensions(gens, n):
    """All extensions of parent (length n) to length n+1 keeping d>=D.
    Returns (list of (r, new gens)), ext counts per r."""
    label, lift, q = quotient_labels(gens, n)
    covered = {label(e) for e in light_errors(n, D - 2)}
    deep = [a for a in range(1 << q) if a not in covered]
    deepset = set(deep)
    out = {0: [list(gens)], 1: [], 2: []}
    top = 2 * n
    for a in deep:
        out[1].append(list(gens) + [lift(a) | (1 << top)])
    ordered_pairs = 0
    for i, a in enumerate(deep):
        for b in deep:
            if b != a and (a ^ b) in deepset:
                ordered_pairs += 1
                if b > a:
                    out[2].append(list(gens) + [lift(a) | (1 << top),
                                                lift(b) | (2 << top)])
    ext = {0: 1, 1: 3 * len(deep), 2: ordered_pairs}
    return out, ext


def to_interleaved(v, n=10):
    """Committed layout (X bits 0..n-1, Z bits n..2n-1) -> interleaved symbols."""
    out = 0
    for i in range(n):
        out |= ((v >> i) & 1) << (2 * i)
        out |= ((v >> (n + i)) & 1) << (2 * i + 1)
    return out


def census(maxn=10, target_k=10):
    cells = {(4, 0): [([], aut_order(graph([], 4)))]}
    assert cells[(4, 0)][0][1] == factorial(4) * 6 ** 4
    for n in range(5, maxn + 1):
        for k in range(max(0, target_k - 2 * (maxn - n)), target_k + 1):
            classes = {}
            lhs = Fraction(0)
            for r in (0, 1, 2):
                for gens, aut in cells.get((n - 1, k - r), []):
                    exts, ext = extensions(gens, n - 1)
                    lhs += Fraction(factorial(n - 1) * 6 ** (n - 1), aut) * ext[r]
                    for g in exts[r]:
                        G = graph(g, n)
                        c = certificate(G)
                        if c not in classes:
                            classes[c] = (g, aut_order(G))
            found = list(classes.values())
            for g, a in found:
                assert (factorial(n) * 6 ** n) % a == 0
                ws = [wt(v, n) for v in span(g) if v]
                assert not ws or min(ws) >= D
                assert len(span(g)) == 1 << k
            rhs = sum(Fraction(factorial(n) * 6 ** n, a) for _, a in found)
            print(f"cell n={n} k={k}: classes={len(found)} "
                  f"labelled mass {lhs} = {rhs}: {lhs == rhs}", flush=True)
            assert lhs == rhs
            cells[(n, k)] = found
    return cells[(maxn, target_k)]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--generators",
                    default=Path(__file__).resolve().parent / "s7_37_classes.json")
    args = ap.parse_args()
    found = census()
    assert len(found) == 37
    committed = [[to_interleaved(v) for v in g]
                 for g in json.loads(Path(args.generators).read_text())]
    mine = {certificate(graph(g, 10)) for g, _ in found}
    theirs = {certificate(graph(g, 10)) for g in committed}
    assert len(theirs) == 37 and mine == theirs
    print("the 37 classes found are the 37 committed representatives")
    auts = [aut_order(graph(g, 10)) for g in committed]
    for i, a in enumerate(auts, 1):
        print(f"class {i:2d}: |Aut| = {a}")
    mass = sum(Fraction(factorial(10) * 6 ** 10, a) for a in auts)
    assert mass == 2333060795842560
    print(f"sum over the 37 classes of 10! 6^10 / |Aut| = {mass}")
    print("ORBIT-STABILIZER MASS CHECK PASS")
