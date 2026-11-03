"""Exact unary/pair filter, with the same semantics as general/xstage.py.

PASS means the necessary conditions survived, not that a prefix map exists.
The CLI distinguishes exhaustive class exclusion from limits and errors. No
running search or submitted manuscript is modified by importing this module.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "general"))
sys.path.insert(0, str(ROOT / "release16" / "search"))
from full import prefix
from sympl import setup
from zstage import leaves
from xstage import solve_affine

D = 5


def _reduce_basis(gens):
    """Independent GF(2) rows, in descending pivot order."""
    piv = {}
    for value in gens:
        value = int(value)
        for bit in sorted(piv, reverse=True):
            if value >> bit & 1:
                value ^= piv[bit]
        if value:
            piv[value.bit_length() - 1] = value
    return tuple((bit, piv[bit]) for bit in sorted(piv, reverse=True))


def projected_affine(p, null, idx, b):
    """Return offset and a basis, without enumerating the projected space."""
    mask = (1 << b) - 1

    def project(value):
        return sum(((value >> (k * b)) & mask) << (t * b)
                   for t, k in enumerate(idx))

    return project(p), _reduce_basis(project(value) for value in null)


def affine_members(values, offset, rows):
    """Exact membership of integer vectors in offset + span(rows)."""
    rem = np.asarray(values, dtype=np.int64) ^ offset
    for bit, row in rows:
        rem = np.where((rem >> bit) & 1, rem ^ row, rem)
    return rem == 0


@lru_cache(maxsize=32)
def _requirements(m, distance, costfn, zbits, xbits):
    """Prefix requirements are independent of the shortened code and Z leaf."""
    a = m - 1
    zwords = [0] * (1 << a)
    for subset in range(1, 1 << a):
        low = subset & -subset
        zwords[subset] = zwords[subset ^ low] ^ zbits[low.bit_length() - 1]

    def requirements(xword):
        return tuple((subset, distance - costfn(xword ^ zword))
                     for subset, zword in enumerate(zwords)
                     if distance - costfn(xword ^ zword) > 0)

    unary = tuple(requirements(xbits[k]) for k in range(a))
    pair = {(k, l): requirements(xbits[k] ^ xbits[l])
            for k in range(a) for l in range(k + 1, a)}
    return unary, pair


def _subset_syndromes(xs):
    syndromes = [0] * (1 << len(xs))
    for subset in range(1, len(syndromes)):
        low = subset & -subset
        syndromes[subset] = syndromes[subset ^ low] ^ int(xs[low.bit_length() - 1])
    return syndromes


def _allowed(dmin, ar, syndromes, requirements):
    allowed = np.ones(len(dmin), dtype=bool)
    for subset, need in requirements:
        allowed &= dmin[ar ^ syndromes[subset]] >= need
        if not allowed.any():
            break
    return allowed


def screen_leaf(xs, dmin, Gam, b, m, costfn, Zbits, Xbits, pairs=True):
    """Return LIN, U, P or PASS, exactly as the original necessary filter."""
    if b > 30:
        raise ValueError("This implementation uses int64 joint vectors; require b <= 30.")
    if len(xs) != m - 1 or len(dmin) != 1 << b:
        raise ValueError("Inconsistent leaf or syndrome dimensions")
    xs = [int(x) for x in xs]
    a = m - 1
    eqs = []
    for j in range(b):
        for l in range(j + 1, b):
            mask = 0
            for i in range(a):
                if xs[i] >> j & 1:
                    mask ^= 1 << (i * b + l)
                if xs[i] >> l & 1:
                    mask ^= 1 << (i * b + j)
            eqs.append((mask, int(Gam[j, l])))
    sol = solve_affine(eqs, a * b)
    if sol is None:
        return "LIN"
    p, null = sol
    dmin = np.asarray(dmin)
    ar = np.arange(1 << b, dtype=np.int64)
    syndromes = _subset_syndromes(xs)
    unary_requirements, pair_requirements = _requirements(
        m, D, costfn, tuple(Zbits), tuple(Xbits))
    domains = []
    for k in range(a):
        allowed = _allowed(dmin, ar, syndromes, unary_requirements[k])
        offset, rows = projected_affine(p, null, [k], b)
        dom = ar[allowed & affine_members(ar, offset, rows)]
        if not len(dom):
            return "U"
        domains.append(dom)
    if not pairs:
        return "PASS"

    # Fixed binary relations on the initial unary domains, built lazily so an
    # early empty pair does not force construction of the remaining pairs.
    # Subsequent pruning only removes rows/columns, so no test is repeated.
    relations = {}
    pair_indices = tuple((k, l) for k in range(a) for l in range(k + 1, a))
    alive = [np.ones(len(dom), dtype=bool) for dom in domains]
    changed = True
    while changed:
        changed = False
        for k, l in pair_indices:
            if (k, l) not in relations:
                dk, dl = domains[k], domains[l]
                allowed = _allowed(dmin, ar, syndromes, pair_requirements[k, l])
                offset, rows = projected_affine(p, null, [k, l], b)
                joint = dk[:, None] | (dl[None, :] << b)
                relations[k, l] = (allowed[dk[:, None] ^ dl[None, :]]
                                   & affine_members(joint, offset, rows))
            relation = relations[k, l]
            active_relation = relation & alive[k][:, None] & alive[l][None, :]
            nk, nl = active_relation.any(axis=1), active_relation.any(axis=0)
            if not nk.any() or not nl.any():
                return "P"
            if not np.array_equal(nk, alive[k]) or not np.array_equal(nl, alive[l]):
                alive[k], alive[l] = nk, nl
                changed = True
    return "PASS"


def run(gens, f, m, pairs=True, *, max_leaves=None, seconds=None,
        progress_every=0, progress=None):
    """Return an explicit completion status and tally. Time limits are soft:
    they are checked between yielded leaves, not inside Z enumeration.
    """
    start = time.monotonic()
    dmin, _, Gam, b = setup(gens, f)
    basis, cost = prefix(m)
    tally = Counter()
    count = 0
    status = "COMPLETE"
    for xs in leaves(dmin, m, d=D):
        if ((max_leaves is not None and count >= max_leaves)
                or (seconds is not None and time.monotonic() - start >= seconds)):
            status = "LIMIT"
            break
        result = screen_leaf(xs, dmin, Gam, b, m, cost, basis[:m - 1], basis[m - 1:], pairs)
        tally[result] += 1
        count += 1
        if progress and progress_every and count % progress_every == 0:
            progress({"status": "RUNNING", "leaves": count, "tally": dict(tally),
                      "seconds": time.monotonic() - start})
    return {"status": status, "excluded": status == "COMPLETE" and not tally["PASS"],
            "leaves": count, "tally": dict(tally), "seconds": time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("m", type=int)
    parser.add_argument("ids", nargs="?", help="Comma-separated class IDs; default all")
    parser.add_argument("--max-leaves", type=int)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()
    if args.n not in (16, 17):
        parser.error("This CLI currently supports the existing n=16,17 distance-5 cells only")
    key = {6: "10,10", 7: "9,8", 8: "8,6", 9: "7,4"}.get(args.m - (args.n - 16))
    if key is None:
        parser.error("No existing cell for these parameters")
    cells = json.loads((ROOT / "release16" / "census" / "cells.json").read_text())[key]
    ids = range(len(cells)) if args.ids is None else [int(i) for i in args.ids.split(",")]
    for i in ids:
        try:
            report = run(cells[i], args.n - args.m, args.m,
                         max_leaves=args.max_leaves, seconds=args.seconds,
                         progress_every=args.progress_every,
                         progress=lambda item: print(json.dumps({"class": i, **item}), flush=True))
            print(json.dumps({"class": i, **report}), flush=True)
        except Exception as error:
            print(json.dumps({"class": i, "status": "ERROR", "excluded": False,
                              "error": f"{type(error).__name__}: {error}"}), flush=True)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
