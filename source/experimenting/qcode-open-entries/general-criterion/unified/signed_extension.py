"""Exact fixed-Z rank-one extension solver, with explicit search-limit status.

Inputs include a verified suffix code via ExtensionModel. Covers 3<=m<=9,
1<=d<=n, b<=14, full-rank Z columns, R=<Z^m>. Does not classify C, localize low-weight
stabilizers, or enumerate all choices of Z columns.
"""
import numpy as np
from model import coordinates, lift_coordinate, reduced_basis
from model import span, weight
from kernel_bound import suffix_stabilizer_bound
from two_site import normal_form, requirements, allowed


class SearchLimit(Exception):
    pass


def full_weight_check(xs, ys, distances, m, d):
    """Check every relevant even-X sector, including four-X at distance five."""
    p = m - 1
    zsyn, ysyn = [0] * (1 << p), [0] * (1 << p)
    for v in range(1, 1 << p):
        low = v & -v
        i = low.bit_length() - 1
        zsyn[v] = zsyn[v ^ low] ^ xs[i]
        ysyn[v] = ysyn[v ^ low] ^ ys[i]
    for x in range(1 << p):
        support = x | ((x.bit_count() & 1) << p)
        for z, need in requirements(m, support, d):
            if x == 0 and z == 0:
                continue
            if distances[ysyn[x] ^ zsyn[z]] < need:
                return False
    return True


def solve_extension(model, xs, *, node_limit=None, minimum_stabilizer_weight=None):
    """Return FOUND/EXCLUDED/LIMIT with exact localized-extension semantics.

    FOUND includes a matrix in model.E coordinates whose full criterion passes.
    EXCLUDED applies only to this C,R and fixed X; it is not a parameter bound.
    Invalid or unsupported data raise ValueError; callers must report ERROR.
    """
    m, d, b = model.m, model.d, model.b
    p = m - 1
    if not (3 <= m <= 9 and 1 <= d <= model.n and b <= 14):
        raise ValueError('Require 3<=m<=9,1<=d<=n,b<=14')
    if model.r != 1 or model.R != [((1 << m) - 1) << m]:
        raise ValueError('Require the normalized local stabilizer R=<Z^m>')
    if len(xs) != p or any(not isinstance(x, int) or not 0 <= x < 1 << b for x in xs):
        raise ValueError('Require m-1 valid integer syndrome columns')
    if node_limit is not None and (not isinstance(node_limit, int) or node_limit < 0):
        raise ValueError('Require a nonnegative integer node limit')
    if minimum_stabilizer_weight is not None and not 1 <= minimum_stabilizer_weight <= m:
        raise ValueError('Minimum stabilizer weight must be between one and m')
    if model.C_distance is not None and model.C_distance < d:
        return dict(status='EXCLUDED', reason='suffix_distance', nodes=0)
    if len(reduced_basis(xs)) != p:
        raise ValueError('Fixed-Z solver requires full column rank')
    nf = normal_form(xs, model.Gamma)
    if any(nf['Gamma'][i] >> j & 1 for i in range(p, b) for j in range(i + 1, b)):
        return dict(status='EXCLUDED', reason='symplectic_Z', nodes=0)
    rank_h = len(reduced_basis([offset >> p for offset in nf['offsets']]))
    suffix_stab_dim = b - p - rank_h
    suffix_distance = max(d, minimum_stabilizer_weight or 1)
    if suffix_stab_dim > suffix_stabilizer_bound(model.f, suffix_distance):
        return dict(status='EXCLUDED', reason='suffix_stabilizer_dimension', nodes=0,
                    suffix_stabilizer_dimension=suffix_stab_dim,
                    minimum_stabilizer_weight_assumption=minimum_stabilizer_weight)
    forward = np.asarray(nf['forward'], dtype=np.int64)
    zsyn = forward[:1 << p]
    pure = requirements(m, 0, d)
    if any(z != 0 and model.distances[zsyn[z]] < need for z, need in pure):
        return dict(status='EXCLUDED', reason='pure_Z_weight', nodes=0)
    distances = np.asarray(model.distances)
    gauges = [(1 << p) - 1 ^ (1 << i) for i in range(p)]
    values = []
    for i in range(p):
        candidates = np.array([u for u in range(1 << p) if not (u >> i & 1)], dtype=np.int64)
        ys = forward[candidates | nf['offsets'][i]]
        domain = candidates[allowed(ys, distances, zsyn,
                                    requirements(m, (1 << i) | (1 << p), d))]
        if not len(domain):
            return dict(status='EXCLUDED', reason='unary_weight', nodes=0)
        assert {int(u) ^ gauges[i] for u in domain} == set(map(int, domain))
        values.append(domain)
    # A global pivot-phase gauge flips every orbit orientation. Fix the first
    # variable's orientation to its smaller representative, retaining all orbits.
    domains = [set(range(len(v))) for v in values]
    domains[0] = {j for j, u in enumerate(values[0]) if u < (int(u) ^ gauges[0])}
    relations = {}

    def relation(i, j):
        if (i, j) not in relations:
            vi, vj = values[i], values[j]
            weight_ok = allowed(np.arange(1 << b), distances, zsyn,
                                requirements(m, (1 << i) | (1 << j), d))
            yi = forward[vi | nf['offsets'][i]]
            yj = forward[vj | nf['offsets'][j]]
            parity = (((vi[:, None] >> j) ^ (vj[None, :] >> i)) & 1)
            relations[i, j] = (weight_ok[yi[:, None] ^ yj[None, :]]
                               & (parity == ((nf['Gamma'][i] >> j) & 1)))
        return relations[i, j]

    def propagate(dom):
        changed = True
        while changed:
            changed = False
            for i in range(p):
                for j in range(i + 1, p):
                    rel = relation(i, j)
                    ni = {u for u in dom[i] if any(rel[u, v] for v in dom[j])}
                    nj = {v for v in dom[j] if any(rel[u, v] for u in dom[i])}
                    if not ni or not nj:
                        return False
                    if ni != dom[i] or nj != dom[j]:
                        dom[i], dom[j] = ni, nj
                        changed = True
        return True

    nodes = rejected_full = 0

    def convert(ys):
        std = ([1 << (m + i) for i in range(p)]
               + [(1 << i) | (1 << p) for i in range(p)])
        return [lift_coordinate(coordinates(e, model.R + std) >> 1, list(xs) + ys)
                for e in model.E]

    def search(dom):
        nonlocal nodes, rejected_full
        nodes += 1
        if node_limit is not None and nodes > node_limit:
            raise SearchLimit
        if not propagate(dom):
            return None
        unresolved = [i for i in range(p) if len(dom[i]) > 1]
        if not unresolved:
            top = [int(values[i][next(iter(dom[i]))]) for i in range(p)]
            ys = [int(forward[u | nf['offsets'][i]]) for i, u in enumerate(top)]
            if full_weight_check(xs, ys, model.distances, m, d):
                if (minimum_stabilizer_weight is not None
                        and minimum_stabilizer_weight > min(m, d)):
                    S, _ = model.construct(convert(ys))
                    if any(weight(v, model.n) < minimum_stabilizer_weight
                           for v in span(S) if v):
                        rejected_full += 1
                        return None
                return top, ys
            rejected_full += 1
            return None
        i = min(unresolved, key=lambda j: (len(dom[j]), j))
        for state in sorted(dom[i]):
            child = [set(v) for v in dom]
            child[i] = {state}
            found = search(child)
            if found is not None:
                return found
        return None

    try:
        found = search(domains)
    except SearchLimit:
        return dict(status='LIMIT', nodes=nodes, node_limit=node_limit,
                    full_weight_rejections=rejected_full)
    if found is None:
        return dict(status='EXCLUDED', reason='signed_CSP_exhausted', nodes=nodes,
                    full_weight_rejections=rejected_full)
    top, ys = found
    columns = convert(ys)
    assert model.criterion(columns)
    return dict(status='FOUND', nodes=nodes, Y=ys, top_columns=top,
                model_columns=columns, full_weight_rejections=rejected_full,
                verified_full_criterion=True,
                minimum_stabilizer_weight_assumption=minimum_stabilizer_weight)
