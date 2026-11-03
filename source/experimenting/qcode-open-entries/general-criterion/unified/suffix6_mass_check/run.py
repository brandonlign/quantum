"""Exact incidence mass from independently checked C5 to all C6 suffixes."""
import json
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from suffix_mass_check.common import canonical, coset_minima
from model import reduced_basis, span, weight
from mixed_frames.stage3.common import atomic_json, sha


def main():
    started = time.monotonic()
    sources = [ROOT / 'opus_parallel/suffix9_dimension5/suffix9_dim5.json',
               ROOT / 'suffix9_classes.json']
    parent = ROOT / 'suffix_mass_check/verification9.json'
    proof = json.loads(parent.read_text())
    assert proof['status'] == 'PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS'
    for name, digest in proof['hashes'].items():
        assert sha(Path(name)) == digest, name
    frozen = {str(p.relative_to(ROOT)): sha(p) for p in sources + [parent]}
    rows5, rows6 = [json.loads(p.read_text())['classes'] for p in sources]
    assert len(rows5) == 32865 and len(rows6) == 56
    command = ['c++', '-O3', '-std=c++17', '-Wall', '-Wextra', str(HERE / 'count.cpp'), '-o', str(HERE / 'count')]
    subprocess.run(command, check=True)
    data = str(len(rows5)) + '\n' + ''.join(str(r['class_index']) + ' ' + ' '.join(map(str, r['representative']['generators'])) + '\n' for r in rows5)
    native = subprocess.run([str(HERE / 'count')], input=data, capture_output=True, text=True, check=True)
    counts = {int(i): (int(n), int(checksum)) for i, n, checksum in [line.split() for line in native.stdout.splitlines()]}
    assert len(counts) == len(rows5)
    # Independent physical-word enumeration; sample includes extendible codes.
    sample = sorted(set(range(8)) | {r['class_index'] for r in rows5 if counts[r['class_index']][0] > 0})
    if len(sample) > 24:
        nonzero = [i for i in sample if i >= 8]
        sample = sorted(set(range(8)) | set(nonzero[::max(1, len(nonzero) // 16)][:16]))
    by_id = {r['class_index']: r for r in rows5}
    for i in sample:
        assert counts[i][0] == int((coset_minima(by_id[i]['representative']['generators'], 9) >= 6).sum()), i
    group = math.factorial(9) * 6 ** 9
    records5, records6 = [], []
    seen5, seen6 = set(), set()
    class_by_certificate, mass5, weighted_hyperplanes = {}, {}, Counter()
    incidence = listed_mass = 0
    for number, r in enumerate(rows5):
        C = r['representative']['generators']
        assert len(reduced_basis(C)) == 5 and min(weight(w, 9) for w in span(C) if w) >= 6
        key, order = canonical(C, 9, need_order=True)
        assert key not in seen5
        seen5.add(key)
        class_by_certificate[key] = r['class_index']
        mass = group // order
        assert mass * order == group
        mass5[r['class_index']] = mass
        n, checksum = counts[r['class_index']]
        incidence += mass * n
        records5.append(dict(cls=r['class_index'], automorphisms=order, orbit_size=mass,
                             extensions=n, quotient_reachability_checksum=checksum))
        if number % 3000 == 0:
            print(json.dumps(dict(stage='dimension5_mass', classes=number, incidence=incidence)), flush=True)
    for r in rows6:
        C = r['representative']['generators']
        assert len(reduced_basis(C)) == 6 and min(weight(w, 9) for w in span(C) if w) >= 6
        key, order = canonical(C, 9, need_order=True)
        assert key not in seen6
        seen6.add(key)
        mass = group // order
        assert mass * order == group
        listed_mass += mass
        records6.append(dict(cls=r['class_index'], automorphisms=order, orbit_size=mass))
        # A second, pointwise check: enumerate every hyperplane of every listed C6.
        for functional in range(1, 64):
            pivot = (functional & -functional).bit_length() - 1
            H = [C[j] ^ (C[pivot] if (functional >> j) & 1 else 0)
                 for j in range(6) if j != pivot]
            assert len(reduced_basis(H)) == 5
            key5 = canonical(H, 9)
            assert key5 in class_by_certificate
            weighted_hyperplanes[class_by_certificate[key5]] += mass
    pointwise = all(weighted_hyperplanes[i] == mass5[i] * n for i, (n, _) in counts.items())
    assert frozen == {str(p.relative_to(ROOT)): sha(p) for p in sources + [parent]}
    match = incidence == 63 * listed_mass and pointwise
    report = dict(status='PASS_RELATIVE_TO_CERTIFIED_C5_CENSUS' if match else 'MASS_MISMATCH_OPEN',
                  group_order=group, dimension5_classes=len(rows5), dimension6_classes=len(rows6),
                  incidence=incidence, listed_dimension6_mass=listed_mass,
                  hyperplanes_per_dimension6_code=63, missing_incidence=incidence - 63 * listed_mass,
                  independent_coset_count_samples=sample, records5=records5, records6=records6,
                  all_32865_pointwise_incidence_counts_match=pointwise,
                  independently_enumerated_hyperplanes=63 * len(rows6),
                  build_command=command, seconds=time.monotonic() - started,
                  hashes=dict(frozen, **{str(p.relative_to(ROOT)): sha(p) for p in
                              [Path(__file__), HERE / 'count.cpp', HERE / 'count', HERE / 'PROOF.md',
                               ROOT / 'suffix_mass_check/common.py', ROOT / 'model.py', ROOT / 'mixed_frames/frames.py']}),
                  scope='C6 classification completeness relative to independently checked C5 census and its recorded release parents; no quantum target exclusion')
    atomic_json(HERE / 'verification.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('records5', 'records6', 'hashes')}), flush=True)
    assert match, 'Listed C6 orbits do not exhaust the incidence mass'


if __name__ == '__main__':
    main()
