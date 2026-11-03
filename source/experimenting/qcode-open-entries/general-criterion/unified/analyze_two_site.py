"""Analyze complete class36 certificates using phase-orbit pair relations."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from model import ExtensionModel
from two_site import normal_form, requirements, allowed

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def analyze():
    C = json.loads((ROOT / 'release16/census/cells.json').read_text())['7,4'][36]
    model = ExtensionModel(16, 5, 5, 9, [511 << 9], C)
    distances = np.asarray(model.distances)
    counts, orbit_sizes, partitions, example = Counter(), Counter(), Counter(), {}
    started = time.monotonic()
    certfile = BASE / 'm9_class36_certificates.jsonl'
    for line in certfile.read_text().splitlines():
        record = json.loads(line)
        if record['verdict'] != 'P':
            continue
        xs, cert = record['xs'], record['certificate']
        nf = normal_form(xs, model.Gamma)
        forward = np.asarray(nf['forward'])
        zsyn = forward[:256]
        domains = cert['initial_domains']
        gauges = [255 ^ (1 << i) for i in range(8)]
        orbit_domains = []
        for i, dom in enumerate(domains):
            assert {u ^ gauges[i] for u in dom} == set(dom)
            reps = np.array(sorted({min(u, u ^ gauges[i]) for u in dom}))
            assert 2 * len(reps) == len(dom)
            orbit_domains.append(reps)
            orbit_sizes[len(reps)] += 1
        colors = [value >> 8 for value in nf['offsets']] + [0]
        partition = tuple(sorted(Counter(colors).values(), reverse=True))
        assert partition[0] >= 3
        partitions[str(partition)] += 1
        relations = {}
        distance_relations = {}

        def pair(i, j):
            if (i, j) in relations:
                return relations[i, j]
            di, dj = orbit_domains[i], orbit_domains[j]
            signs = (((di[:, None] >> j) ^ (dj[None, :] >> i)
                      ^ (nf['Gamma'][i] >> j)) & 1)
            weight_ok = allowed(np.arange(1024), distances, zsyn,
                                requirements(9, (1 << i) | (1 << j), 5))
            yi = forward[di | nf['offsets'][i]]
            yj0 = forward[dj | nf['offsets'][j]]
            yj1 = forward[(dj ^ gauges[j]) | nf['offsets'][j]]
            good0 = weight_ok[yi[:, None] ^ yj0[None, :]]
            good1 = weight_ok[yi[:, None] ^ yj1[None, :]]
            relation = np.where(signs == 0, good0, good1)
            # Simultaneously flipping both orientations preserves every weight.
            yi1 = forward[(di ^ gauges[i]) | nf['offsets'][i]]
            assert np.array_equal(good0, weight_ok[yi1[:, None] ^ yj1[None, :]])
            assert np.array_equal(good1, weight_ok[yi1[:, None] ^ yj0[None, :]])
            relations[i, j] = relation
            distance_relations[i, j] = good0 | good1
            return relation

        direct = None
        for i in range(8):
            for j in range(i + 1, 8):
                if not pair(i, j).any():
                    direct = (i, j)
                    break
            if direct is not None:
                break
        if direct is not None:
            counts['initial_empty_pair'] += 1
            kind = ('compatibility_essential' if distance_relations[direct].any()
                    else 'distance_empty_within_fixed_cosets')
            counts[kind] += 1
            counts['initial_empty_monochromatic_pivot_triple'] += int(
                colors[direct[0]] == colors[direct[1]] == 0)
            example.setdefault(kind, dict(xs=xs, empty_pair=list(direct),
                                         orbit_sizes=[len(v) for v in orbit_domains],
                                         colors=colors))
        else:
            counts['requires_propagation'] += 1
            example.setdefault('requires_propagation', dict(xs=xs, colors=colors,
                                                           orbit_sizes=[len(v) for v in orbit_domains]))
        alive = [np.ones(len(dom), dtype=bool) for dom in orbit_domains]
        empty = False
        changed = True
        while changed and not empty:
            changed = False
            for i in range(8):
                for j in range(i + 1, 8):
                    active = pair(i, j) & alive[i][:, None] & alive[j][None, :]
                    ni, nj = active.any(axis=1), active.any(axis=0)
                    if not ni.any() or not nj.any():
                        empty = True
                        break
                    if not np.array_equal(ni, alive[i]) or not np.array_equal(nj, alive[j]):
                        alive[i], alive[j] = ni, nj
                        changed = True
                if empty:
                    break
        assert empty, 'Orbit quotient and original AC verdict disagree'
        counts['orbit_filter_agreements'] += 1
        if counts['orbit_filter_agreements'] % 250 == 0:
            print(json.dumps(dict(status='RUNNING', checked=counts['orbit_filter_agreements'],
                                  seconds=time.monotonic() - started)), flush=True)
    assert counts['orbit_filter_agreements'] == 1101
    report = dict(status='PASS', class_index=36, counts=dict(counts),
                  orbit_domain_sizes=dict(sorted(orbit_sizes.items())),
                  color_partitions=dict(partitions), examples=example,
                  seconds=time.monotonic() - started,
                  certificate_sha256=hashlib.sha256(certfile.read_bytes()).hexdigest(),
                  analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Analysis of all1101 pair-rejected leaves of the completed class; '
                        'phase-orbit relations preserve AC rejection; no general color bound proved')
    (BASE / 'm9_class36_structure.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('examples', 'color_partitions')}))


if __name__ == '__main__':
    analyze()
