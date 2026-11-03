"""Independent complete-leaf comparisons for Gamma and unary-X pruning.

The oracle evaluates physical prefix weights over R=<Z^m> directly. It does
not copy the C unary cost formula or the incremental domain/rank algorithm.
"""
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from model import ExtensionModel, embed, span, weight
from run_benchmark import z_condition
from signed_extension import solve_extension


def run(binary, m, d, distances, gamma, flags):
    payload = f'{m} {len(distances)} {d}\n' + ' '.join(map(str, distances)) + '\n'
    if '--gamma' in flags:
        payload += ' '.join(map(str, gamma)) + '\n'
    p = subprocess.run([str(HERE / binary), *flags], input=payload,
                       text=True, capture_output=True, timeout=30, check=True)
    stats = json.loads(p.stderr)
    assert stats['status'] == 'COMPLETE', stats
    leaves = [tuple(map(int, line.split()[1:])) for line in p.stdout.splitlines()]
    assert len(leaves) == stats['leaves']
    assert len(leaves) == len(set(leaves))
    return stats, leaves


def unary_ok(xs, distances, m, d):
    """Exhaust every suffix image for each physical X_i X_pivot sector."""
    translates = []
    for mask in range(1 << (m - 1)):
        syndrome = 0
        for j, x in enumerate(xs):
            if mask >> j & 1:
                syndrome ^= x
        translates.append(syndrome)
    for i in range(m - 1):
        physical_x = (1 << i) | (1 << (m - 1))
        costs = [min(weight(physical_x | (z << m), m),
                     weight(physical_x | ((z ^ ((1 << m) - 1)) << m), m))
                 for z in range(1 << (m - 1))]
        if not any(all(cost + distances[y ^ syndrome] >= d
                       for cost, syndrome in zip(costs, translates))
                   for y in range(len(distances))):
            return False
    return True


def verify():
    rng = random.Random(1736)
    cases = []
    for case in range(80):
        m = 3 + case % 2
        b = 3 + case % 2
        distances = [0] + [rng.randint(1, 4) for _ in range((1 << b) - 1)]
        gamma = [0] * b
        for i in range(b):
            for j in range(i):
                if rng.randrange(2):
                    gamma[i] ^= 1 << j
                    gamma[j] ^= 1 << i
        s0, plain = run('z_enum8192', m, 3, distances, gamma, [])
        expected = sorted(xs for xs in plain
                          if z_condition(xs, gamma, b) and unary_ok(xs, distances, m, 3))
        s1, actual = run('z_enum_xdom', m, 3, distances, gamma, ['--gamma', '--xdom'])
        assert sorted(actual) == expected, (case, expected, actual)
        cases.append(dict(case=case, plain=len(plain), survivors=len(actual),
                          plain_nodes=s0['nodes'], pruned_nodes=s1['nodes']))

    # C=0 on two suffix sites, R=<ZZZ>: accepted map reconstructs a real
    # five-qubit stabilizer state with minimum stabilizer weight >=3.
    model = ExtensionModel(5, 0, 3, 3, [7 << 3], [])
    s0, plain = run('z_enum8192', 3, 3, model.distances, model.Gamma, [])
    expected = sorted(xs for xs in plain if z_condition(xs, model.Gamma, model.b)
                      and unary_ok(xs, model.distances, 3, 3))
    s1, actual = run('z_enum_xdom', 3, 3, model.distances, model.Gamma,
                     ['--gamma', '--xdom'])
    assert actual and sorted(actual) == expected
    positive = None
    for xs in actual:
        answer = solve_extension(model, list(xs), minimum_stabilizer_weight=3)
        if answer['status'] == 'FOUND':
            S, D = model.construct(answer['model_columns'])
            minimum = min(weight(w, 5) for w in span(S) if w)
            assert len(S) == 5 and minimum >= 3 and set(span(S)) == set(span(D))
            positive = dict(xs=xs, generators=S, minimum_weight=minimum)
            break
    assert positive, 'pruning must retain a genuine positive code'
    # Zero node budgets cannot become complete exclusions.
    payload = '3 16 3\n' + ' '.join(map(str, model.distances)) + '\n'
    payload += ' '.join(map(str, model.Gamma)) + '\n'
    limited = subprocess.run([str(HERE/'z_enum_xdom'), '--gamma', '--xdom',
                              '--max-nodes', '0'], input=payload,
                             text=True, capture_output=True, check=True)
    assert json.loads(limited.stderr)['status'] == 'LIMIT'
    result = dict(status='PASS', cases=cases,
                  nonempty_cases=sum(c['survivors'] > 0 for c in cases),
                  positive_control=positive,
                  positive_plain_leaves=len(plain), positive_pruned_leaves=len(actual),
                  scope='Complete random-table comparisons and a physical positive state; not a target exclusion',
                  sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in ['z_enum_xdom.c', 'z_enum_xdom', 'verify_xdom_prune.py']})
    (HERE/'xdom_prune_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases', 'sha256')}))


if __name__ == '__main__':
    verify()
