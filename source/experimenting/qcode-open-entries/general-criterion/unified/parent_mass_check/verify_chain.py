"""Bind independent parent, lengthening, lifting and extension certificates."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from mixed_frames.stage3.common import atomic_json, sha


def main():
    names = ['parent_mass_check/verification.json', 'suffix_mass_check/verification8.json',
             'suffix_mass_check/verification9.json', 'suffix6_mass_check/verification.json']
    reports = []
    for name in names:
        path = ROOT / name
        report = json.loads(path.read_text())
        for source, digest in report['hashes'].items():
            assert sha(ROOT / source) == digest, (name, source)
        reports.append(report)
    parents, eight, nine, six = reports
    assert parents['status'] == 'PASS_COMPLETE_7_SITE_D5_PARENTS_FROM_TYPES'
    assert eight['status'] == nine['status'] == 'PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS'
    assert six['status'] == 'PASS_RELATIVE_TO_CERTIFIED_C5_CENSUS'
    assert parents['dimension3_type_enumeration'] == eight['parent_dimension3']
    assert all(r['every_parent_pointwise_incidence_matches'] for r in parents['transitions'])
    for k in (3, 4, 5):
        p = ROOT / f'opus_parallel/suffix9_dimension5/cells_d5_7_{k}.json'
        assert parents['hashes'][str(p.relative_to(ROOT))] == eight['hashes'][str(p)]
    assert nine['parent_verification_sha256'] == sha(ROOT / names[1])
    p8 = ROOT / 'opus_parallel/suffix9_dimension5/cell_8_5.json'
    assert eight['hashes'][str(p8)] == nine['hashes'][str(p8)]
    p9 = ROOT / 'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    assert nine['hashes'][str(p9)] == six['hashes'][str(p9.relative_to(ROOT))]
    assert six['hashes']['suffix_mass_check/verification9.json'] == sha(ROOT / names[2])
    assert eight['children'] == nine['parent_classes'] == 26299
    assert nine['source_classes'] == six['dimension5_classes'] == 32865
    assert six['dimension6_classes'] == 56 and six['missing_incidence'] == 0
    assert six['all_32865_pointwise_incidence_counts_match']
    report = dict(status='PASS_EXHAUSTIVE_SUFFIX_CLASSIFICATION_CHAIN',
                  cells=[dict(sites=7, dimension=3, distance=5, classes=32),
                         dict(sites=7, dimension=4, distance=5, classes=43),
                         dict(sites=7, dimension=5, distance=5, classes=1),
                         dict(sites=8, dimension=5, distance=5, classes=26299),
                         dict(sites=9, dimension=5, distance=6, classes=32865),
                         dict(sites=9, dimension=6, distance=6, classes=56)],
                  inherited_release_parent_completeness_discharged=True,
                  dependencies=names,
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in
                          [ROOT / name for name in names] + [Path(__file__), HERE / 'CHAIN.md']},
                  scope='Exhaustive additive-code classification coverage, using independent coordinate types, pointwise incidences, lengthening mass and complete puncture lifts. This supplies classification premises only; [[17,3,6]] extension coverage and final theorem review remain open.')
    atomic_json(HERE / 'chain.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'hashes'}), flush=True)


if __name__ == '__main__':
    main()
