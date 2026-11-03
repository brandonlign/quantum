"""Exact fixed-Z completion with free diagonal and all complementary sectors.

No hull-zero assumption, coefficient-space labeling, or diagonal gauge.
COMPLETE means exhaustion for this fixed labeled Z only. FOUND and LIMIT are
not exhaustion. Every reported solution passes model.criterion and a physical
centralizer-span check independent of the cost/sector tables.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from model import (ExtensionModel, coordinates, dual, embed, lift_coordinate,
                   nullspace, reduced_basis, symp)
from .frames import mixed_frame, xor_table

_HERE = Path(__file__).resolve().parent
_HASHES = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [Path(__file__), _HERE/'frames.py', _HERE.parent/'model.py']}


def normal_form(xs, Gamma):
    """Syndrome congruence, following two_site.normal_form without its imports."""
    p, b = len(xs), len(Gamma)
    if len(reduced_basis(xs)) != p or any(x < 0 or x >= 1 << b for x in xs):
        raise ValueError('Fixed Z columns must be independent syndrome labels')
    T = list(xs)
    for j in range(b):
        if len(reduced_basis(T + [1 << j])) > len(T):
            T.append(1 << j)
    invcols = [coordinates(1 << j, T) for j in range(b)]
    invrows = [sum(((c >> i) & 1) << j for j, c in enumerate(invcols)) for i in range(b)]
    transformed = []
    for u in invrows:
        transformed.append(sum(((u & lift_coordinate(v, Gamma)).bit_count() & 1) << j
                               for j, v in enumerate(invrows)))
    offsets = [sum(((transformed[p+j] >> i) & 1) << (p+j)
                   for j in range(b-p)) for i in range(p)]
    return dict(p=p, b=b, T=T, inverse_rows=invrows, Gamma=transformed,
                offsets=offsets, forward=xor_table(T))


def physical_check(model, frame, columns):
    """Construct S as the physical adjoint graph, and enumerate its dual.

    Does not use prefix costs, suffix distances, sector predicates or the
    normal form. Integer masks still use the repository Pauli convention.
    """
    p = frame.p
    S = [embed(r, 0, model.m, model.f) for r in frame.R]
    for j, h in enumerate(model.P):
        pre = 0
        for i in range(p):
            if columns[p+i] >> j & 1:
                pre ^= frame.E[i]
            if columns[i] >> j & 1:
                pre ^= frame.E[p+i]
        S.append(embed(pre, h, model.m, model.f))
    S = reduced_basis(S)
    assert len(S) == model.n-model.k
    assert all(symp(u, v, model.n) == 0 for u in S for v in S)
    D = dual(S, model.n)
    assert len(D) == model.n+model.k
    model_columns = frame.model_columns(model, columns)
    oldS, oldD = model.construct(model_columns)
    assert reduced_basis(oldS) == S and reduced_basis(oldD) == reduced_basis(D)
    # Shortening via coefficient nullspaces (no span-weight oracle).
    def shorten(words, part):
        from model import split
        projections = [split(w, model.m, model.f)[1-part] for w in words]
        width = 2*(model.f if part == 0 else model.m)
        rows = [sum(((w >> j) & 1) << i for i, w in enumerate(projections))
                for j in range(width)]
        return reduced_basis([split(lift_coordinate(c, words), model.m, model.f)[part]
                              for c in nullspace(rows, len(words))])
    assert shorten(S, 0) == frame.R
    assert shorten(D, 1) == model.C
    Rwords = set(xor_table([embed(r, 0, model.m, model.f) for r in frame.R]))
    mask = (1 << model.n)-1
    minimum_outside = model.n+1
    # Gray traversal holds one physical word, rather than a 2^(n+k) list.
    word = 0
    pure_distance = model.n+1
    for a in range(1 << len(D)):
        if a:
            low = a & -a
            word ^= D[low.bit_length()-1]
        wt = ((word & mask) | (word >> model.n)).bit_count()
        if word:
            pure_distance = min(pure_distance, wt)
        if word not in Rwords:
            minimum_outside = min(minimum_outside, wt)
    assert minimum_outside >= model.d
    return dict(S=S, D=D, centralizer_words_checked=1 << len(D),
                minimum_outside_R=minimum_outside, pure_distance=pure_distance)


def validate_solution(model, frame, columns):
    assert len(columns) == 2*frame.p
    images = xor_table(columns)
    assert all(frame.costs[e]+model.distances[images[e]] >= model.d
               for e in range(1, len(images)))
    model_columns = frame.model_columns(model, columns)
    assert model.criterion(model_columns)
    physical = physical_check(model, frame, columns)
    return dict(adapted_columns=list(columns), model_columns=model_columns,
                kernel=model.kernel_rank_identity(model_columns), **physical)


class _Limit(Exception):
    pass


def solve(model, frame, xs, *, find_all=False, node_limit=None,
          seconds=None, solution_limit=None):
    """Enumerate exact Y solutions for an independent fixed labeled X.

    Limits cover preprocessing as well as DFS. A solution cap is LIMIT in
    find_all mode, even if it happens to equal the eventual solution count.
    """
    start = time.monotonic()
    if len(xs) != frame.p:
        raise ValueError('Wrong number of fixed Z columns')
    if model.m != frame.m or model.R != frame.R:
        raise ValueError('Frame/model mismatch')
    for limit in (node_limit, seconds, solution_limit):
        if limit is not None and limit < 0:
            raise ValueError('Limits must be nonnegative')
    nodes, leaves, solutions = 0, 0, []
    reason = None
    p, Q = frame.p, 1 << frame.p

    def check_limit():
        if seconds is not None and time.monotonic()-start >= seconds:
            raise _Limit('time')

    def result(status, complete):
        return dict(status=status, complete=complete,
                    fixed_z_excluded=complete and not solutions,
                    reason=reason, X=list(xs), solutions=solutions,
                    nodes=nodes, leaves=leaves, seconds=time.monotonic()-start,
                    limits=dict(nodes=node_limit, seconds=seconds, solutions=solution_limit),
                    context=dict(n=model.n, k=model.k, d=model.d, m=model.m,
                                 frame=frame.name, R=frame.R, E=frame.E, C=model.C,
                                 P=model.P, hull=model.hull,
                                 label_convention='s_j=<P_j,w>; e=z|(x<<p)'),
                    implementation_sha256=_HASHES,
                    scope='This fixed labeled Z only; no class or frame exclusion.')

    try:
        check_limit()
        if node_limit == 0:
            raise _Limit('nodes')
        if solution_limit == 0:
            raise _Limit('solutions')
        nf = normal_form(xs, model.Gamma)
        G = nf['Gamma']
        if any(G[i] >> p for i in range(p, model.b)):
            reason = 'nonzero_annihilator_Gram'
            return result('COMPLETE', True)
        if model.C_distance is not None and model.C_distance < model.d:
            reason = 'suffix_distance'
            return result('COMPLETE', True)
        zsyn = xor_table(xs)
        if any(frame.costs[z]+model.distances[zsyn[z]] < model.d for z in range(1, Q)):
            reason = 'pure_Z_weight'
            return result('COMPLETE', True)
        offset = xor_table(nf['offsets'])
        # allowed[S] is a bitset of all low-coordinate XOR values, diagonal
        # included. No sharing between sectors with different physical costs.
        allowed = [0] * Q
        for sector in range(1, Q):
            check_limit()
            requirements = [(z, model.d-frame.costs[z | (sector << p)])
                            for z in range(Q)
                            if frame.costs[z | (sector << p)] < model.d]
            for u in range(Q):
                beta = nf['forward'][u ^ offset[sector]]
                if all(model.distances[beta ^ zsyn[z]] >= need for z, need in requirements):
                    allowed[sector] |= 1 << u
            if not allowed[sector]:
                reason = 'empty_sector_'+str(sector)
                return result('COMPLETE', True)
        domains = [allowed[1 << i] for i in range(p)]
        values = [None] * p

        def translate(bits, shift):
            out = 0
            while bits:
                low = bits & -bits
                out |= 1 << ((low.bit_length()-1) ^ shift)
                bits ^= low
            return out

        def dfs(domains, assigned):
            nonlocal nodes, leaves
            check_limit()
            if node_limit is not None and nodes >= node_limit:
                raise _Limit('nodes')
            nodes += 1
            # Every sector whose last unassigned variable is v constrains its
            # domain by translating allowed[S] by the assigned XOR.
            domains = domains[:]
            for sector in range(1, Q):
                remaining = sector & ~assigned
                total = 0
                done = sector & assigned
                while done:
                    bit = done & -done
                    total ^= values[bit.bit_length()-1]
                    done ^= bit
                if not remaining:
                    if not (allowed[sector] >> total & 1):
                        return False
                elif remaining & (remaining-1) == 0:
                    i = remaining.bit_length()-1
                    domains[i] &= translate(allowed[sector], total)
                    if not domains[i]:
                        return False
            if assigned == Q-1:
                leaves += 1
                assert all((((values[i] >> j) ^ (values[j] >> i)) & 1)
                           == ((G[i] >> j) & 1) for i in range(p) for j in range(i+1, p))
                ys = [nf['forward'][values[i] ^ nf['offsets'][i]] for i in range(p)]
                solution = validate_solution(model, frame, list(xs)+ys)
                solution['U'] = values[:]
                solutions.append(solution)
                check_limit()
                if not find_all:
                    return True
                if solution_limit is not None and len(solutions) >= solution_limit:
                    raise _Limit('solutions')
                return False
            i = min((j for j in range(p) if not assigned >> j & 1),
                    key=lambda j: domains[j].bit_count())
            options = domains[i]
            while options:
                bit = options & -options
                u = bit.bit_length()-1
                options ^= bit
                values[i] = u
                new = domains[:]
                feasible = True
                # All off-diagonal U_ji+U_ij=G_ij equations; U_ii untouched.
                for j in range(p):
                    if j == i or assigned >> j & 1:
                        continue
                    want = ((G[i] >> j) ^ (u >> j)) & 1
                    new[j] &= sum(1 << v for v in range(Q) if (v >> i & 1) == want)
                    if not new[j]:
                        feasible = False
                        break
                if feasible and dfs(new, assigned | (1 << i)):
                    values[i] = None
                    return True
                values[i] = None
            return False

        if dfs(domains, 0):
            return result('FOUND', False)
        return result('COMPLETE', True)
    except _Limit as exc:
        reason = str(exc)
        return result('LIMIT', False)
    except KeyboardInterrupt:
        reason = 'interrupt'
        return result('LIMIT', False)


def main():
    from .suffix_tables import SOURCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frame', choices=['668', '688', '778'], required=True)
    parser.add_argument('--class-index', type=int, required=True)
    parser.add_argument('--z', type=json.loads, required=True)
    parser.add_argument('--nodes', type=int)
    parser.add_argument('--seconds', type=float)
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()
    c = next(c for c in json.loads(SOURCE.read_text())['classes']
             if c['class_index'] == args.class_index)
    frame = mixed_frame(args.frame)
    model = ExtensionModel(17, 3, 6, 8, frame.R, c['representative']['generators'])
    r = solve(model, frame, args.z, find_all=args.all,
              node_limit=args.nodes, seconds=args.seconds)
    r.update(class_index=args.class_index, mu=6,
             source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest())
    print(json.dumps(r))


if __name__ == '__main__':
    main()
