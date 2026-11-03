"""Bounded, explicit-status canonical Z + unary/pair benchmark.

Does not replace or terminate the original workers. COMPLETE/EXCLUDED requires
exhaustion; a single surviving leaf establishes only failure of this filter.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from z_enumerator import ZEnumerator, EnumerationLimit
from z_c_backend import CEnumerator
from x_filter_fast import screen_leaf
from model import ExtensionModel, nullspace, weight

ROOT = Path(__file__).resolve().parents[2]


def z_condition(xs, Gamma, b):
    B = nullspace(xs, b)
    for c in B:
        image = sum(((row & c).bit_count() & 1) << i
                    for i, row in enumerate(Gamma))
        if any((v & image).bit_count() & 1 for v in B):
            return False
    return True


def run_class(C, n, m, *, seconds=None, max_nodes=None, stage='x', progress_every=1000,
              backend='python'):
    started = time.monotonic()
    p = ((1 << m) - 1) << m
    model = ExtensionModel(n, n - 11, 5, m, [p], C)
    enum_type = CEnumerator if backend == 'c' else ZEnumerator
    enumerator = enum_type(model.distances, m, seconds=seconds, max_nodes=max_nodes)
    Gamma = np.array([[(row >> j) & 1 for j in range(model.b)]
                      for row in model.Gamma], dtype=np.uint8)
    mask = (1 << m) - 1
    physical_cost = lambda v: min(weight(v, m), weight(v ^ p, m))
    Zbits = [1 << (m + i) for i in range(m - 1)]
    Xbits = [(1 << i) | (1 << (m - 1)) for i in range(m - 1)]
    tally = Counter()
    status = 'EXCLUDED'
    leaves = enumerator.leaves()
    try:
        for xs in leaves:
            if not z_condition(xs, model.Gamma, model.b):
                result = 'Z'
            elif stage == 'z':
                result = 'PASS'
            else:
                result = screen_leaf(xs, model.distances, Gamma, model.b, m,
                                     physical_cost, Zbits, Xbits)
                assert result != 'LIN', 'Z condition and affine compatibility disagree'
            tally[result] += 1
            if result == 'PASS':
                status = 'SURVIVES'
                break
            if progress_every and enumerator.leaf_count % progress_every == 0:
                print(json.dumps({'status': 'RUNNING', 'leaves': enumerator.leaf_count,
                                  'nodes': enumerator.nodes if backend == 'python' else None,
                                  'tally': dict(tally),
                                  'seconds': time.monotonic() - started}), flush=True)
    except EnumerationLimit as error:
        status = 'LIMIT'
        reason = str(error)
    finally:
        leaves.close()
    report = {'status': status, 'excluded': status == 'EXCLUDED',
              'completed_exclusion': status == 'EXCLUDED',
              'canonical_Z': True, 'stage': stage, 'backend': backend,
              'leaves': enumerator.leaf_count,
              'nodes': (None if backend == 'c' and enumerator.summary is None
                        else enumerator.nodes), 'tally': dict(tally),
              'seconds': time.monotonic() - started}
    if status == 'LIMIT':
        report['limit_reason'] = reason
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, choices=[16, 17])
    parser.add_argument('m', type=int)
    parser.add_argument('ids', nargs='?')
    parser.add_argument('--stage', choices=['z', 'x'], default='x')
    parser.add_argument('--backend', choices=['python', 'c'], default='python')
    parser.add_argument('--seconds', type=float)
    parser.add_argument('--max-nodes', type=int)
    parser.add_argument('--progress-every', type=int, default=1000)
    args = parser.parse_args()
    key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}.get(args.m - (args.n - 16))
    if key is None:
        parser.error('No existing classified cell for these parameters')
    data = ROOT / 'release16/census/cells.json'
    cells = json.loads(data.read_text())[key]
    ids = range(len(cells)) if args.ids is None else [int(i) for i in args.ids.split(',')]
    files = [data, Path(__file__), Path(__file__).with_name('model.py'),
             Path(__file__).with_name('z_enumerator.py'),
             Path(__file__).with_name('x_filter_fast.py'), ROOT / 'general/xstage.py']
    if args.backend == 'c':
        files += [Path(__file__).with_name('z_c_backend.py'),
                  Path(__file__).with_name('z_enumerator.c'),
                  Path(__file__).with_name('z_enumerator_c')]
    print(json.dumps({'status': 'INPUTS', 'n': args.n, 'k': args.n - 11, 'd': 5,
                      'm': args.m, 'cell': key,
                      'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in files}}), flush=True)
    for i in ids:
        print(json.dumps({'class': i, 'status': 'STARTED'}), flush=True)
        try:
            report = run_class(cells[i], args.n, args.m, seconds=args.seconds,
                               max_nodes=args.max_nodes, stage=args.stage,
                               progress_every=args.progress_every, backend=args.backend)
        except Exception as error:
            print(json.dumps({'class': i, 'status': 'ERROR', 'excluded': False,
                              'error': f'{type(error).__name__}: {error}'}), flush=True)
            return 1
        print(json.dumps({'class': i, **report}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
