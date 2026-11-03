"""Nonempty controls and exhaustive solution/leaf-set comparisons.

python3 -B -m mixed_frames.verify --all-controls
All checks always run; the flag documents the intended invocation.
Results are written only under mixed_frames/.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import random
import time

from model import (ExtensionModel, dual, embed, lift_coordinate, quotient_basis,
                   nullspace, pullback_gram, recover_extension, reduced_basis, span, symp, weight)
from .frames import FRAMES, Frame, adapted_frame, mixed_frame, xor_table
from .fixed_z import normal_form, solve, validate_solution
from .enumerate_z import enumerate_z, restricted_rank
from .suffix_tables import verify_physical, build

HERE = Path(__file__).resolve().parent


def physical_span(words):
    """Independent un-reduced span implementation for physical oracles."""
    result = {0}
    for w in words:
        result |= {a ^ w for a in result}
    return result


def physical_weight(w, n):
    return sum(bool(((w >> i) | (w >> (n+i))) & 1) for i in range(n))


def physical_reps(model):
    reps = {}
    for w in range(1 << (2*model.f)):
        # Independent per-site commutation spelling, no model.syndrome.
        syndrome = 0
        for j, h in enumerate(model.P):
            bit = sum((((h >> i) & 1)*((w >> (model.f+i)) & 1)
                       + ((w >> i) & 1)*((h >> (model.f+i)) & 1))
                      for i in range(model.f)) & 1
            syndrome |= bit << j
        reps.setdefault(syndrome, w)
    assert len(reps) == 1 << model.b
    return reps


def brute_profile(model, frame, xs):
    """All physical Y tuples; no U, normal form, compatibility or cost table.

    The adjoint graph is formed in physical Pauli space, its commutators are
    checked, and its centralizer is independently enumerated. This tests the
    complete Y set, including every free diagonal choice.
    """
    profile = {}
    p = frame.p
    Rglobal = {embed(r, 0, model.m, model.f) for r in physical_span(frame.R)}
    count = 0
    for ys in itertools.product(range(1 << model.b), repeat=p):
        count += 1
        columns = list(xs) + list(ys)
        S = [embed(r, 0, model.m, model.f) for r in frame.R]
        for j, h in enumerate(model.P):
            v = 0
            for i in range(p):
                if ys[i] >> j & 1:
                    v ^= frame.E[i]
                if xs[i] >> j & 1:
                    v ^= frame.E[p+i]
            S.append(embed(v, h, model.m, model.f))
        if any(symp(u, v, model.n) for i, u in enumerate(S) for v in S[:i]):
            continue
        assert len(physical_span(S)) == 1 << (model.n-model.k)
        # Separate physical coset graph for D, using representatives obtained
        # by enumerating all physical suffix words.
        reps = physical_reps(model)
        D = [embed(r, 0, model.m, model.f) for r in frame.R]
        D += [embed(0, c, model.m, model.f) for c in model.C]
        D += [embed(e, reps[s], model.m, model.f) for e, s in zip(frame.E, columns)]
        words = physical_span(D)
        assert len(words) == 1 << (model.n+model.k)
        assert all(symp(w, s, model.n) == 0 for w in D for s in S)
        profile[ys] = min(physical_weight(w, model.n) for w in words if w not in Rglobal)
    return profile, count


def independent_xs(p, b):
    return [xs for xs in itertools.product(range(1, 1 << b), repeat=p)
            if len(physical_span(xs)) == 1 << p]


def brute_z_leaves(model, frame, rank_filter, sector_filter):
    """Exhaustive labeled X oracle using physical prefix/suffix fibers.

    Rank feasibility is tested by existence of *some* compatible physical
    adjoint graph over every Y tuple, not by a restricted-rank budget.
    """
    p = frame.p
    reps = physical_reps(model)
    codewords = physical_span(model.C)
    delta = [min(physical_weight(reps[s] ^ c, model.f) for c in codewords)
             for s in range(1 << model.b)]
    Rwords = physical_span(frame.R)
    costs = [min(physical_weight(w ^ r, frame.m) for r in Rwords)
             for w in xor_table(frame.E)]
    leaves = []
    for xs in independent_xs(p, model.b):
        images = xor_table(xs)
        if any(costs[z]+delta[images[z]] < model.d for z in range(1, 1 << p)):
            continue
        if rank_filter:
            profile, _ = brute_profile(model, frame, xs)
            if not profile:
                continue
        if sector_filter:
            if any(not any(all(costs[z | (s << p)]+delta[beta ^ images[z]] >= model.d
                               for z in range(1 << p)) for beta in range(1 << model.b))
                   for s in range(1, 1 << p) if s.bit_count() <= 2):
                continue
        leaves.append(list(xs))
    return leaves


def compare_small(model, frame, thresholds, label, samples=None):
    rows, matrices = [], 0
    controls = 0
    for xs in samples if samples is not None else independent_xs(frame.p, model.b):
        profile, count = brute_profile(model, frame, xs)
        matrices += count
        for d in thresholds:
            model.d = d
            expected = {ys for ys, distance in profile.items() if distance >= d}
            result = solve(model, frame, list(xs), find_all=True)
            assert result['status'] == 'COMPLETE' and result['complete'], result
            got = {tuple(s['adapted_columns'][frame.p:]) for s in result['solutions']}
            assert len(got) == len(result['solutions'])
            assert got == expected, (label, xs, d, expected, got)
            assert result['fixed_z_excluded'] == (not expected)
            controls += bool(got)
            rows.append(dict(X=list(xs), d=d, solutions=len(got)))
    return dict(label=label, cases=len(rows), nonempty_cases=controls,
                physical_Y_tuples=matrices, rows=rows)


def normalize_on_word(S, n, word):
    """Physical local symplectic maps send the chosen word to Z^m."""
    support = [i for i in range(n) if ((word & ((1 << n)-1)) | (word >> n)) >> i & 1]
    order = support + [i for i in range(n) if i not in support]
    out = []
    for w in S:
        x, z = 0, 0
        for j, i in enumerate(order):
            a, b = (w >> i) & 1, (w >> (n+i)) & 1
            if j < len(support):
                px, pz = (word >> i) & 1, (word >> (n+i)) & 1
                u = (a*pz) ^ (b*px)
                # Complement X for a Z/Y word, Z for an X word.
                qx, qz = (1, 0) if pz else (0, 1)
                ra, rb = a ^ (u*qx), b ^ (u*qz)
                assert (ra, rb) in ((0, 0), (px, pz))
                a, b = u, int(bool(ra or rb))
            x |= a << j
            z |= b << j
        out.append(x | (z << n))
    return out, len(support)


def control_from_stabilizer(S, n, k, d, label, word=None):
    if word is None:
        word = min((w for w in physical_span(S) if w), key=lambda w: (physical_weight(w, n), w))
    normalized, m = normalize_on_word(S, n, word)
    model, known = recover_extension(normalized, n, k, d, m)
    frame = adapted_frame(m, model.R, label)
    columns = frame.adapted_columns(model, known)
    original = validate_solution(model, frame, columns)
    result = solve(model, frame, columns[:frame.p], node_limit=200000, seconds=60)
    assert result['status'] == 'FOUND', result
    assert result['solutions'][0]['pure_distance'] >= d
    nf = normal_form(columns[:frame.p], model.Gamma)
    U = [sum(((row & y).bit_count() & 1) << j for j, row in enumerate(nf['inverse_rows']))
         & ((1 << frame.p)-1) for y in columns[frame.p:]]
    return dict(label=label, n=n, k=k, d=d, m=m, hull=model.hull,
                known=original, known_U=U,
                known_diagonal=[(u >> i) & 1 for i, u in enumerate(U)], solver=result)


def positive_controls():
    def pauli(s):
        return sum(int(c in 'XY') << i for i, c in enumerate(s)) | (
            sum(int(c in 'YZ') << (len(s)+i) for i, c in enumerate(s)))
    fiveS = [pauli(s) for s in ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']]
    five = control_from_stabilizer(fiveS, 5, 1, 3, '[[5,1,3]]', word=fiveS[0])
    assert any(five['known_diagonal'])
    model = ExtensionModel(7, 1, 2, 4, [240, 15], [24, 3])
    frame = Frame(4, model.R, model.E, model.cost, 'hull-two supplied basis')
    known = validate_solution(model, frame, [1, 2, 0, 8])
    assert known['kernel'] == (1, 0, 2) and model.hull == 2
    assert known['pure_distance'] == 2
    assert reduced_basis(known['S']) == reduced_basis([8197, 8202, 1328, 832, 1920, 6144])
    r = solve(model, frame, [1, 2], find_all=True)
    assert r['status'] == 'COMPLETE'
    assert any(s['adapted_columns'] == [1, 2, 0, 8] for s in r['solutions'])
    controls = [five, dict(label='n=7 hull-two', known=known, solver=r,
                           basis_note='Supplied model.E; physical pure-Z X is dependent at d=2.')]
    for name in FRAMES:
        frame = mixed_frame(name)
        model = ExtensionModel(14, 0, 2, 8, frame.R, [])
        columns = [model.syndrome(1 << (6+i)) for i in range(6)] + [model.syndrome(1 << i) for i in range(6)]
        known = validate_solution(model, frame, columns)
        assert known['pure_distance'] == 2
        r = solve(model, frame, columns[:6])
        assert r['status'] == 'FOUND' and r['solutions'][0]['pure_distance'] == 2
        controls.append(dict(label='n=14 mixed graph '+name, known=known, solver=r))
    # Hull-six physical fixture at measured distance two, not target distance six.
    frame = mixed_frame('688')
    C = [((1 << i) | (1 << 8)) << 9 for i in range(6)]
    model = ExtensionModel(17, 3, 2, 8, frame.R, C)
    suffixE = quotient_basis(model.C, model.f)
    xs = [model.syndrome(w) for w in suffixE[:3]]
    for j in range(model.b):
        if len(xs) == 6:
            break
        if len(reduced_basis(xs + [1 << j])) > len(xs):
            xs.append(1 << j)
    r = solve(model, frame, xs)
    assert r['status'] == 'FOUND' and model.hull == 6
    assert r['solutions'][0]['pure_distance'] == 2
    assert r['solutions'][0]['kernel'] == (6, 6, 6)
    controls.append(dict(label='physical hull-six fixture, d=2', solver=r))
    # Stored sixteen-qubit distance-six control, without invoking its writer.
    source = HERE.parent / 'positive_control_16_1_6.txt'
    S = []
    for line in source.read_text().splitlines():
        x, z = line.split('|')
        S.append(sum(int(c) << i for i, c in enumerate(x.split())) | (
            sum(int(c) << (16+i) for i, c in enumerate(z.split()))))
    controls.append(control_from_stabilizer(S, 16, 1, 6, 'stored [[16,1,6]]'))
    return controls


def limits_and_interrupts(model, frame, xs):
    import mixed_frames.fixed_z as fixed
    results = []
    for kw in [dict(node_limit=0), dict(node_limit=1), dict(seconds=0),
               dict(find_all=True, solution_limit=0), dict(find_all=True, solution_limit=1)]:
        r = solve(model, frame, xs, **kw)
        assert r['status'] == 'LIMIT' and not r['complete'] and not r['fixed_z_excluded']
        results.append(dict(kind='fixed_Z', limits=kw, reason=r['reason']))
    original = fixed.validate_solution
    def interrupt(*args):
        raise KeyboardInterrupt
    fixed.validate_solution = interrupt
    try:
        r = solve(model, frame, xs)
        assert r['status'] == 'LIMIT' and r['reason'] == 'interrupt'
    finally:
        fixed.validate_solution = original
    results.append(dict(kind='fixed_Z', reason='interrupt'))
    for kw in [dict(node_limit=0), dict(seconds=0), dict(leaf_limit=0), dict(leaf_limit=1)]:
        r = enumerate_z(model, frame, **kw)
        assert r['status'] == 'LIMIT' and not r['complete']
        results.append(dict(kind='Z_enumerator', limits=kw, reason=r['reason']))
    r = enumerate_z(model, frame, on_leaf=interrupt)
    assert r['status'] == 'LIMIT' and r['reason'] == 'interrupt'
    results.append(dict(kind='Z_enumerator', reason='interrupt'))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-controls', action='store_true')
    parser.parse_args()
    started = time.monotonic()
    report = dict(status='PASS', controls=positive_controls(), exhaustive=[], Z_comparisons=[])
    print('PASS: nonempty physical controls (5, 7, three mixed 14, hull-six, stored 16)', flush=True)
    # Full-plane n=6 control: all 210 labeled independent X maps and all 256 Y
    # tuples for each. Includes the nonempty [[6,0,4]] code.
    model = ExtensionModel(6, 0, 1, 4, [15, 240], [])
    frame = adapted_frame(4, model.R, 'full-plane n6')
    row = compare_small(model, frame, [1, 2, 3, 4], 'n6 full-plane all independent X')
    assert any(r['d'] == 4 and r['solutions'] for r in row['rows'])
    report['exhaustive'].append(row)
    print('PASS: n=6 exhaustive fixed-Z/Y comparison, including [[6,0,4]]', flush=True)
    model.d = 3
    for rank_filter, sector_filter in [(False, False), (True, False), (False, True), (True, True)]:
        expected = brute_z_leaves(model, frame, rank_filter, sector_filter)
        r = enumerate_z(model, frame, use_rank=rank_filter, use_domains=sector_filter)
        assert r['status'] == 'COMPLETE' and r['leaves'] == expected
        report['Z_comparisons'].append(dict(label='n6 d3', rank=rank_filter,
                                            sectors=sector_filter, leaves=r['leaf_count']))
    all_z = enumerate_z(model, frame)
    report['limits'] = limits_and_interrupts(model, frame, all_z['leaves'][0])
    # Degenerate b=4: the same complete tuple comparison, not hull-zero labels.
    model = ExtensionModel(7, 1, 1, 4, [240, 15], [24, 3])
    frame = Frame(4, model.R, model.E, model.cost, 'hull-two')
    report['exhaustive'].append(compare_small(model, frame, [1, 2, 3], 'n7 hull-two all independent X'))
    print('PASS: n=7 hull-two exhaustive fixed-Z/Y comparison and limit controls', flush=True)
    # Every length-two binary dimension-two suffix subspace, all hull regimes.
    subspaces = {tuple(reduced_basis([a, b])) for a in range(1, 16) for b in range(a+1, 16)
                 if len(reduced_basis([a, b])) == 2}
    assert len(subspaces) == 35
    for C in sorted(subspaces):
        model = ExtensionModel(4, 1, 1, 2, [12], list(C))
        frame = adapted_frame(2, model.R)
        report['exhaustive'].append(compare_small(model, frame, [1, 2, 3], 'n4 C='+str(C)))
        for d in (2, 3):
            model.d = d
            expected = brute_z_leaves(model, frame, True, True)
            r = enumerate_z(model, frame)
            assert r['status'] == 'COMPLETE' and r['leaves'] == expected
            report['Z_comparisons'].append(dict(label='n4 C='+str(C), d=d, leaves=r['leaf_count']))
    # Three- and four-variable fixtures exercise sectors larger than pairs.
    # Every Y tuple is checked for each listed X, rather than sampling Y.
    model = ExtensionModel(6, 1, 1, 4, [240], [])
    frame = adapted_frame(4, model.R, 'p3')
    report['exhaustive'].append(compare_small(model, frame, [1, 2, 3, 4],
        'p3 all Y, four fixed X', samples=[(1, 2, 4), (3, 5, 9), (7, 11, 13), (15, 14, 12)]))
    model = ExtensionModel(7, 2, 1, 5, [31 << 5], [])
    frame = adapted_frame(5, model.R, 'p4')
    report['exhaustive'].append(compare_small(model, frame, [1, 2, 3],
        'p4 all Y, one fixed X', samples=[(1, 2, 4, 8)]))
    print('PASS: exhaustive Y for three/four-variable higher-sector fixtures', flush=True)
    # Check normal-form basis changes against the original compatibility
    # matrix for arbitrary U, retaining diagonal, on every target suffix.
    from .suffix_tables import SOURCE
    rng = random.Random(20261013)
    report['normal_form_checks'] = 0
    report['physical_partial_rank_checks'] = 0
    for c in json.loads(SOURCE.read_text())['classes']:
        t = build(c['representative']['generators'])
        xs = []
        while len(xs) < 6:
            candidate = rng.randrange(1, 4096)
            if len(reduced_basis(xs+[candidate])) > len(xs):
                xs.append(candidate)
        for j in range(7):
            ann = nullspace(xs[:j], 12)
            assert len(ann) == 12-j
            assert all((a & x).bit_count() % 2 == 0 for a in ann for x in xs[:j])
            physical_ann = [lift_coordinate(a, t['P']) for a in ann]
            physical_gram = [sum(symp(a, b, 9) << k for k, b in enumerate(physical_ann))
                             for a in physical_ann]
            assert restricted_rank(xs[:j], t['Gamma']) == len(reduced_basis(physical_gram))
            report['physical_partial_rank_checks'] += 1
        # The forward map and the inverse coefficient rows must be inverse.
        nf = normal_form(xs, t['Gamma'])
        assert all(sum(((row & w).bit_count() & 1) << j for j, row in enumerate(nf['inverse_rows'])) == e
                   for e, w in enumerate(nf['forward']))
        for _ in range(2):
            # Independently choose all 21 free bits and fill the other triangle.
            U = [0]*6
            for i in range(6):
                U[i] |= rng.randrange(2) << i
                for j in range(i+1, 6):
                    bit = rng.randrange(2)
                    U[i] |= bit << j
                    U[j] |= (bit ^ ((nf['Gamma'][i] >> j) & 1)) << i
            ys = [nf['forward'][u ^ off] for u, off in zip(U, nf['offsets'])]
            physical_form = pullback_gram(xs+ys, 12)
            transformed = [sum(((a & lift_coordinate(b, physical_form)).bit_count() & 1) << j
                               for j, b in enumerate(nf['inverse_rows'])) for a in nf['inverse_rows']]
            expected = nf['Gamma'][:6] + [g & 63 for g in nf['Gamma'][6:]]
            assert transformed == expected
            report['normal_form_checks'] += 1
    report['seconds'] = time.monotonic()-started
    report['fixed_Z_cases'] = sum(r['cases'] for r in report['exhaustive'])
    report['physical_Y_tuples'] = sum(r['physical_Y_tuples'] for r in report['exhaustive'])
    report['sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in list(HERE.glob('*.py')) + [HERE.parent/'model.py',
                                                          HERE.parent/'positive_control_16_1_6.txt']}
    report['scope'] = ('Exact fixed-Z and labeled necessary-Z small-case solution sets; nonempty physical controls. '
                       'No target distance-six extension exclusion or complete target search.')
    (HERE/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('controls', 'exhaustive', 'sha256', 'Z_comparisons', 'limits')}))


if __name__ == '__main__':
    main()
