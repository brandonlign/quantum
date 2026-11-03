"""Independent unrestricted local lift replay and physical capacity checks."""
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from pynauty import certificate
from model import coordinates, dual, gram, reduced_basis, span, symp, weight
from suffix_mass_check.common import star_graph, coset_minima
from mixed_frames.stage3.common import atomic_json, sha


def distinguished(H, p):
    words = span(H)
    W = len(words)
    graph = star_graph(H, 8)
    mark = words.index(p)
    graph.set_vertex_coloring([{mark}, set(range(W)) - {mark},
                               set(range(W, W + 24)), set(range(W + 24, W + 32))])
    return certificate(graph)


def main():
    started = time.monotonic()
    parent_proof = ROOT / 'parent_mass_check/chain.json'
    chain = json.loads(parent_proof.read_text())
    assert chain['status'] == 'PASS_EXHAUSTIVE_SUFFIX_CLASSIFICATION_CHAIN'
    for name, digest in chain['hashes'].items():
        assert sha(ROOT / name) == digest, name
    plan_path = ROOT / 'kernel_case_plan/plan.json'
    profiles_path = ROOT / 'kernel_capacity/profiles.json'
    plan = json.loads(plan_path.read_text())
    profiles = json.loads(profiles_path.read_text())['records']
    actual, expected = defaultdict(set), defaultdict(set)
    p = 255 << 8
    valid_h = 0
    for h in range(1 << 16):
        if h >= (h ^ p) or h in (0, p) or symp(h, p, 8):
            continue
        if min(weight(h, 8), weight(h ^ p, 8)) < 6:
            continue
        actual[1, 0].add(distinguished([p, h], p))
        valid_h += 1
    assert valid_h == 7072 and len(actual[1, 0]) == 4
    replays = []
    parents_paths = []
    for dim in (3, 4):
        source = ROOT / f'opus_parallel/suffix9_dimension5/cells_d5_7_{dim}.json'
        parents_paths.append(source)
        parents = json.loads(source.read_text())[f'7,{dim}']['xz']
        valid = distinguished_pairs = 0
        for C in parents:
            assert len(reduced_basis(C)) == dim
            assert min(weight(w, 7) for w in span(C) if w) >= 5
            for a in range(1 << dim):
                for b in range(1 << dim):
                    # Independent direct seven-to-eight physical embedding.
                    H = [(v & 127) | ((v >> 7) << 8) |
                         (((a >> j) & 1) << 7) | (((b >> j) & 1) << 15)
                         for j, v in enumerate(C)]
                    words = span(H)
                    if min(weight(w, 8) for w in words if w) < 6:
                        continue
                    valid += 1
                    t = len(reduced_basis(gram(H, 8)))
                    for marked in words:
                        if weight(marked, 8) == 8 and all(symp(marked, w, 8) == 0 for w in H):
                            assert t in (0, 2)
                            actual[dim - 1, t].add(distinguished(H, marked))
                            distinguished_pairs += 1
        replays.append(dict(q=dim - 1, parents=len(parents), ordered_pairs=len(parents) * (1 << (2 * dim)),
                            minimum6_lifts=valid, marked_full_radical_words=distinguished_pairs))
    for prefix in plan['prefix']:
        key = distinguished(prefix['H'], prefix['p'])
        assert key.hex() == prefix['distinguished_certificate']
        assert key not in expected[prefix['q'], prefix['t']]
        expected[prefix['q'], prefix['t']].add(key)
    assert dict(actual) == dict(expected), 'Missing or extra distinguished local kernel class'
    capacity_checks = []
    for i, r in enumerate(profiles):
        H, p = r['H'], r['p']
        assert len(reduced_basis(H)) == r['q'] + 1 and p in span(H)
        assert weight(p, 8) == 8 and all(symp(p, h, 8) == 0 for h in H)
        assert min(weight(w, 8) for w in span(H) if w) >= 6
        assert len(reduced_basis(gram(H, 8))) == r['t']
        assert distinguished(H, p) in actual[r['q'], r['t']]
        quotient = dual(H, 8)
        p_coordinates = coordinates(p, quotient)
        distances = coset_minima(H, 8)
        retained = [int(d) for label, d in enumerate(distances)
                    if (label & p_coordinates).bit_count() % 2 == 0]
        assert len(retained) == 1 << r['quotient_bits'] == 1 << (14 - r['q'])
        histogram = Counter(retained)
        capacities = [sum(d >= 6 - j for d in retained if d) for j in range(6)]
        assert histogram == {int(k): v for k, v in r['histogram'].items()}, i
        assert capacities == r['caps'], i
        capacity_checks.append(dict(profile=i, q=r['q'], t=r['t'], cosets=len(retained), caps=capacities))
    paths = [Path(__file__), HERE / 'PROOF.md', parent_proof, plan_path, profiles_path,
             ROOT / 'model.py', ROOT / 'suffix_mass_check/common.py', ROOT / 'mixed_frames/frames.py'] + parents_paths
    report = dict(status='PASS_COMPLETE_LOCAL_KERNEL_CLASSES_AND_CAPACITIES',
                  q1_direct_normalized_subspaces=valid_h, ordered_lift_replays=replays,
                  distinguished_classes={f'{q},{t}': len(keys) for (q, t), keys in sorted(actual.items())},
                  capacity_profiles=len(capacity_checks), capacity_checks=capacity_checks,
                  seconds=time.monotonic() - started,
                  hashes={str(path.relative_to(ROOT)): sha(path) for path in paths},
                  scope='Complete distinguished local H,p classes and independent physical capacity histograms. No quantum extension or target exclusion; full search coverage remains required.')
    atomic_json(HERE / 'verification.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('hashes', 'capacity_checks')}), flush=True)


if __name__ == '__main__':
    main()
