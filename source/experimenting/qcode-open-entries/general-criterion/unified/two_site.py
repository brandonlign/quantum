"""Independent fixed-Z normal-form filter and replayable rejection certificates.

Uses a syndrome basis change and scalar pair equations, not solve_affine or
projected_affine. Pair consistency is necessary only; PASS is not existence.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from model import ExtensionModel, coordinates, reduced_basis
from z_c_backend import CEnumerator

ROOT = Path(__file__).resolve().parents[2]


def normal_form(xs, Gamma):
    a, b = len(xs), len(Gamma)
    if len(reduced_basis(xs)) != a:
        raise ValueError('The fixed-Z normal form requires independent columns')
    T = list(xs)
    for j in range(b):
        if len(reduced_basis(T + [1 << j])) > len(T):
            T.append(1 << j)
    if len(T) != b:
        raise ValueError('Failed to complete a syndrome basis')
    inverse_columns = [coordinates(1 << j, T) for j in range(b)]
    inverse_rows = [sum(((c >> i) & 1) << j for j, c in enumerate(inverse_columns))
                    for i in range(b)]
    images = [sum(((row & w).bit_count() & 1) << j
                  for j, row in enumerate(Gamma)) for w in inverse_rows]
    transformed = [sum(((w & images[j]).bit_count() & 1) << j
                       for j in range(b)) for w in inverse_rows]
    forward = [0] * (1 << b)
    for v in range(1, 1 << b):
        low = v & -v
        forward[v] = forward[v ^ low] ^ T[low.bit_length() - 1]
    offsets = [sum(((transformed[a + j] >> i) & 1) << (a + j)
                   for j in range(b - a)) for i in range(a)]
    return dict(a=a, b=b, T=T, inverse_rows=inverse_rows,
                Gamma=transformed, offsets=offsets, forward=forward)


@lru_cache(maxsize=None)
def requirements(m, support, d, physical=False):
    result = []
    q = support.bit_count()
    for z in range(1 << (m - 1)):
        if physical:
            all_sites = (1 << m) - 1
            cost = min((support | z).bit_count(),
                       (support | (z ^ all_sites)).bit_count())
        else:
            t = (z & ~support).bit_count()
            cost = q + min(t, m - q - t)
        if d > cost:
            result.append((z, d - cost))
    return tuple(sorted(result, key=lambda item: (-item[1], item[0])))


def allowed(values, distances, z_syndromes, req):
    ok = np.ones(len(values), dtype=bool)
    for z, need in req:
        ok &= distances[values ^ z_syndromes[z]] >= need
        if not ok.any():
            break
    return ok


def screen(xs, distances, Gamma, m, d=5):
    nf = normal_form(xs, Gamma)
    a, b = nf['a'], nf['b']
    transformed = nf['Gamma']
    for i in range(a, b):
        for j in range(i + 1, b):
            if transformed[i] >> j & 1:
                return dict(verdict='Z', certificate=dict(type='Z', i=i, j=j))
    arr = np.asarray(distances)
    forward = np.asarray(nf['forward'], dtype=np.int64)
    zsyn = forward[:1 << a]
    domains = []
    for i in range(a):
        # Diagonal gauge: B_ii=0. Physical weights and all pair tests are invariant.
        us = np.array([u for u in range(1 << a) if not (u >> i & 1)], dtype=np.int64)
        ys = forward[us | nf['offsets'][i]]
        req = requirements(m, (1 << i) | (1 << a), d)
        domain = us[allowed(ys, arr, zsyn, req)]
        if not len(domain):
            return dict(verdict='U', certificate=dict(type='U', variable=i))
        domains.append(domain)
    initial = [list(map(int, values)) for values in domains]
    alive = [np.ones(len(values), dtype=bool) for values in domains]
    relations, trace = {}, []
    changed = True
    while changed:
        changed = False
        for i in range(a):
            for j in range(i + 1, a):
                if (i, j) not in relations:
                    di, dj = domains[i], domains[j]
                    yi = forward[di | nf['offsets'][i]]
                    yj = forward[dj | nf['offsets'][j]]
                    weight_ok = allowed(np.arange(1 << b), arr, zsyn,
                                        requirements(m, (1 << i) | (1 << j), d))
                    parity = (((di[:, None] >> j) ^ (dj[None, :] >> i)) & 1)
                    relations[i, j] = (weight_ok[yi[:, None] ^ yj[None, :]]
                                       & (parity == ((transformed[i] >> j) & 1)))
                active = relations[i, j] & alive[i][:, None] & alive[j][None, :]
                ni, nj = active.any(axis=1), active.any(axis=0)
                if not ni.any() or not nj.any():
                    cert = dict(type='P', initial_domains=initial, trace=trace,
                                empty_pair=[i, j])
                    return dict(verdict='P', certificate=cert)
                if not np.array_equal(ni, alive[i]) or not np.array_equal(nj, alive[j]):
                    trace.append(dict(pair=[i, j],
                                      removed_i=list(map(int, domains[i][alive[i] & ~ni])),
                                      removed_j=list(map(int, domains[j][alive[j] & ~nj]))))
                    alive[i], alive[j] = ni, nj
                    changed = True
    return dict(verdict='PASS', certificate=None)


def verify_certificate(xs, distances, Gamma, m, certificate, d=5):
    """Pure-Python physical-cost replay; does not call screen or its relations.

    Reconstructs the maximal unary domains, checks each justified deletion, and
    verifies the final empty relation. Shares the audited basis transformation.
    """
    nf = normal_form(xs, Gamma)
    a, b = nf['a'], nf['b']
    if certificate['type'] == 'Z':
        c, e = [nf['inverse_rows'][certificate[key]] for key in ('i', 'j')]
        assert all(not (c & x).bit_count() % 2 and not (e & x).bit_count() % 2 for x in xs)
        image = sum(((row & e).bit_count() & 1) << j for j, row in enumerate(Gamma))
        assert (c & image).bit_count() % 2
        return True
    forward, offsets = nf['forward'], nf['offsets']
    zsyn = forward[:1 << a]

    @lru_cache(maxsize=None)
    def admissible(support, syndrome):
        return all(distances[syndrome ^ zsyn[z]] >= need
                   for z, need in requirements(m, support, d, physical=True))

    domains = []
    for i in range(a):
        domain = {u for u in range(1 << a) if not (u >> i & 1)
                  and admissible((1 << i) | (1 << a), forward[u | offsets[i]])}
        domains.append(domain)
    if certificate['type'] == 'U':
        assert not domains[certificate['variable']]
        return True
    assert certificate['type'] == 'P'
    assert domains == list(map(set, certificate['initial_domains']))

    def relation(i, j, u, v):
        if (((u >> j) ^ (v >> i)) & 1) != ((nf['Gamma'][i] >> j) & 1):
            return False
        syndrome = forward[u | offsets[i]] ^ forward[v | offsets[j]]
        return admissible((1 << i) | (1 << j), syndrome)

    for step in certificate['trace']:
        i, j = step['pair']
        remove_i, remove_j = set(step['removed_i']), set(step['removed_j'])
        assert remove_i <= domains[i] and remove_j <= domains[j]
        assert all(not any(relation(i, j, u, v) for v in domains[j]) for u in remove_i)
        assert all(not any(relation(i, j, u, v) for u in domains[i]) for v in remove_j)
        domains[i] -= remove_i
        domains[j] -= remove_j
    i, j = certificate['empty_pair']
    assert not any(relation(i, j, u, v) for u in domains[i] for v in domains[j])
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--class-index', type=int, default=36)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--certificates', type=Path, required=True)
    args = parser.parse_args()
    data = ROOT / 'release16/census/cells.json'
    C = json.loads(data.read_text())['7,4'][args.class_index]
    model = ExtensionModel(16, 5, 5, 9, [511 << 9], C)
    from x_filter_fast import screen_leaf
    Gam = np.array([[(row >> j) & 1 for j in range(model.b)]
                    for row in model.Gamma], dtype=np.uint8)
    Zbits = [1 << (9 + i) for i in range(8)]
    Xbits = [(1 << i) | (1 << 8) for i in range(8)]
    from model import weight
    cost = lambda v: min(weight(v, 9), weight(v ^ (511 << 9), 9))
    enumerator = CEnumerator(model.distances, 9)
    tally, trace_sizes, pairs, example = Counter(), Counter(), Counter(), {}
    started = time.monotonic()
    with args.certificates.open('w') as records:
        for xs in enumerator.leaves():
            result = screen(xs, model.distances, model.Gamma, 9)
            old = screen_leaf(list(xs), model.distances, Gam, model.b, 9, cost, Zbits, Xbits)
            assert result['verdict'] == ('Z' if old == 'LIN' else old)
            cert = result['certificate']
            if cert is None:
                raise AssertionError('Surviving leaf: cannot claim exclusion')
            assert verify_certificate(xs, model.distances, model.Gamma, 9, cert)
            tally[result['verdict']] += 1
            if result['verdict'] == 'P':
                trace_sizes[len(cert['trace'])] += 1
                pairs[tuple(cert['empty_pair'])] += 1
            record = dict(xs=list(xs), **result)
            records.write(json.dumps(record) + '\n')
            example.setdefault(result['verdict'], record)
            if enumerator.leaf_count % 500 == 0:
                print(json.dumps(dict(status='RUNNING', leaves=enumerator.leaf_count,
                                      tally=dict(tally), seconds=time.monotonic() - started)), flush=True)
    assert enumerator.summary['status'] == 'COMPLETE'
    files = [data, Path(__file__), Path(__file__).with_name('model.py'),
             Path(__file__).with_name('z_c_backend.py'),
             Path(__file__).with_name('z_enumerator.c'),
             Path(__file__).with_name('z_enumerator_c'),
             Path(__file__).with_name('x_filter_fast.py'), ROOT / 'general/xstage.py']
    report = dict(status='EXCLUDED', completed_exclusion=True, class_index=args.class_index,
                  leaves=enumerator.leaf_count, nodes=enumerator.nodes, tally=dict(tally),
                  pure_certificate_replays=enumerator.leaf_count,
                  old_filter_agreements=enumerator.leaf_count,
                  pair_trace_lengths=dict(trace_sizes),
                  empty_pairs={str(k): v for k, v in pairs.items()}, examples=example,
                  seconds=time.monotonic() - started,
                  sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in files},
                  certificate_sha256=hashlib.sha256(args.certificates.read_bytes()).hexdigest(),
                  scope='Full canonical class enumeration shared with the original benchmark; '
                        'independent normal-form filter, no old affine solver; '
                        'pure physical-cost certificate replay; not novelty or other-class exclusion')
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('examples', 'sha256')}), flush=True)


if __name__ == '__main__':
    main()
