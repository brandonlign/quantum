"""Reproduce the independent reviewer's deterministic broader finite check.

This compares algebraic conditions with physical construction using model.py;
it is not a separate implementation of the full large-case search.
"""
import argparse
import itertools
import json
from pathlib import Path
import random
from model import ExtensionModel, reduced_basis, symp


def verify():
    rng = random.Random(8675309)
    checks = accepted = 0
    cases = [(3, 1, 2, 1), (4, 1, 3, 2), (5, 1, 3, 2),
             (5, 2, 3, 2), (4, 1, 2, 0), (4, 1, 1, 1)]
    counts = []
    for n, k, m, r in cases:
        before = checks
        f = n - m
        cdim = n + k - 2 * m + r
        for trial in range(12):
            R = []
            while len(R) < r:
                v = rng.randrange(1, 1 << (2 * m))
                if (all(not symp(v, u, m) for u in R)
                        and len(reduced_basis(R + [v])) > len(R)):
                    R = reduced_basis(R + [v])
            C = []
            while len(C) < cdim:
                C = reduced_basis(C + [rng.randrange(1, 1 << (2 * f))])
            for d in range(1, 5):
                model = ExtensionModel(n, k, d, m, R, C)
                for M in itertools.product(range(1 << model.b), repeat=model.a):
                    exact = model.criterion(M)
                    assert exact == model.direct_property(M)
                    model.kernel_rank_identity(M)
                    checks += 1
                    accepted += exact
        counts.append(dict(n=n, k=k, m=m, r=r, instances=checks - before))
    assert checks == 198000 and accepted == 7525
    return dict(status='PASS', instances=checks, accepted=accepted,
                seed=8675309, cases=counts,
                scope='12 sampled data pairs per case (repeats possible), all maps, '
                      'd=1..4; shared model implementation; no novelty or large-search audit')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = verify()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
