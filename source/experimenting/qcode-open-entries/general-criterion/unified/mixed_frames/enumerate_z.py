"""Complete labeled independent-Z reference enumerator, with safe prunes.

No quotient by prefix/suffix symmetry. All labels in [0,2^b) are considered.
Optional complementary singleton/pair domains are weight-only relaxations.
COMPLETE concerns necessary Z leaves, not exact Y completions or exclusions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from model import ExtensionModel, lift_coordinate, nullspace, reduced_basis
from .frames import mixed_frame
from .suffix_tables import SOURCE


def restricted_rank(xs, Gamma):
    ann = nullspace(xs, len(Gamma))
    rows = [sum(((a & lift_coordinate(b, Gamma)).bit_count() & 1) << j
                for j, b in enumerate(ann)) for a in ann]
    return len(reduced_basis(rows))


def enumerate_z(model, frame, *, node_limit=None, seconds=None, leaf_limit=None,
                use_rank=True, use_domains=True, on_leaf=None, store_leaves=True):
    if max(frame.cost[0]) >= model.d:
        raise ValueError('Completeness requires injectivity forced by pure-Z costs')
    if frame.m != model.m or frame.R != model.R:
        raise ValueError('Frame/model mismatch')
    for limit in (node_limit, seconds, leaf_limit):
        if limit is not None and limit < 0:
            raise ValueError('Limits must be nonnegative')
    start = time.monotonic()
    p, Q, N = frame.p, 1 << frame.p, 1 << model.b
    sectors = [s for s in range(1, Q) if s.bit_count() <= 2]
    needs = {s: [(z, model.d-frame.costs[z | (s << p)])
                  for z in range(Q)] for s in sectors}
    leaves, count, nodes = [], 0, 0
    prunes = dict(dependent=0, pure_z=0, rank=0, sector=0)
    reason = None

    class Limit(Exception):
        pass

    def limits():
        if seconds is not None and time.monotonic()-start >= seconds:
            raise Limit('time')

    def dfs(xs, images, domains):
        nonlocal nodes, count
        limits()
        if node_limit is not None and nodes >= node_limit:
            raise Limit('nodes')
        nodes += 1
        j = len(xs)
        if j == p:
            count += 1
            if store_leaves:
                leaves.append(xs[:])
            if on_leaf is not None:
                on_leaf(xs[:])
            if leaf_limit is not None and count >= leaf_limit:
                raise Limit('leaves')
            return
        # No image sorting, circuit condition, or special exclusion of zero
        # complementary syndromes. Independence alone excludes zero Z images.
        for beta in range(N):
            limits()
            if beta in images:
                prunes['dependent'] += 1
                continue
            new_images = [a ^ beta for a in images]
            if any(frame.costs[(1 << j) | z]+model.distances[a] < model.d
                   for z, a in enumerate(new_images)):
                prunes['pure_z'] += 1
                continue
            new_xs = xs + [beta]
            if use_rank and restricted_rank(new_xs, model.Gamma) > 2*(p-j-1):
                prunes['rank'] += 1
                continue
            new_domains = domains
            if use_domains:
                new_domains = {}
                bad = False
                # Retain beta=0 in every complementary domain. Newly assigned
                # pure-Z subsets use their actual labels, never cardinality.
                for sector in sectors:
                    survivors = []
                    for candidate in domains[sector]:
                        limits()
                        if all(model.distances[candidate ^ a] >= needs[sector][(1 << j) | z][1]
                               for z, a in enumerate(new_images)):
                            survivors.append(candidate)
                    if not survivors:
                        bad = True
                        break
                    new_domains[sector] = survivors
                if bad:
                    prunes['sector'] += 1
                    continue
            dfs(new_xs, images + new_images, new_domains)

    try:
        limits()
        if node_limit == 0:
            raise Limit('nodes')
        if leaf_limit == 0:
            raise Limit('leaves')
        initial = {s: [a for a in range(N) if model.distances[a] >= needs[s][0][1]]
                   for s in sectors} if use_domains else {}
        if use_domains and any(not d for d in initial.values()):
            prunes['sector'] += 1
        else:
            dfs([], [0], initial)
        status = 'COMPLETE'
    except Limit as exc:
        status, reason = 'LIMIT', str(exc)
    except KeyboardInterrupt:
        status, reason = 'LIMIT', 'interrupt'
    return dict(status=status, complete=status == 'COMPLETE', reason=reason,
                nodes=nodes, leaf_count=count, leaves=leaves, prunes=prunes,
                seconds=time.monotonic()-start,
                limits=dict(nodes=node_limit, seconds=seconds, leaves=leaf_limit),
                prunes_enabled=dict(rank=use_rank, complementary_domains=use_domains),
                scope='Labeled Z necessary-filter leaves only; Y completion not run; no exclusion claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frame', choices=['668', '688', '778'], required=True)
    parser.add_argument('--class-index', type=int, required=True)
    parser.add_argument('--nodes', type=int)
    parser.add_argument('--seconds', type=float)
    parser.add_argument('--leaves', type=int)
    parser.add_argument('--no-domains', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    if args.output and not args.output.resolve().is_relative_to(here):
        parser.error('Output must be inside mixed_frames/')
    frame = mixed_frame(args.frame)
    c = next(c for c in json.loads(SOURCE.read_text())['classes']
             if c['class_index'] == args.class_index)
    model = ExtensionModel(17, 3, 6, 8, frame.R, c['representative']['generators'])
    r = enumerate_z(model, frame, node_limit=args.nodes, seconds=args.seconds,
                    leaf_limit=args.leaves, use_domains=not args.no_domains)
    r.update(frame=args.frame, class_index=args.class_index, d=6, mu=6,
             E=frame.E, R=frame.R, C=model.C, P=model.P, hull=model.hull,
             label_convention='s_j=<P_j,w>; e=z|(x<<6)',
             source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             implementation_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in [Path(__file__), here/'frames.py', here.parent/'model.py']})
    if args.output:
        args.output.write_text(json.dumps(r, indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k not in ('leaves', 'implementation_sha256')}))


if __name__ == '__main__':
    main()
