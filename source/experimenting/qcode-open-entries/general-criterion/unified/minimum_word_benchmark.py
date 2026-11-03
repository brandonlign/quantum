"""Bounded all43 diagnostic for the minimum-stabilizer suffix obstruction."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
from model import ExtensionModel, reduced_basis
from two_site import normal_form
from z_c_backend import CEnumerator
from z_enumerator import EnumerationLimit

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def main():
    data = ROOT / 'release16/census/cells.json'
    cells = json.loads(data.read_text())['7,4']
    reports = []
    started = time.monotonic()
    for i, C in enumerate(cells):
        model = ExtensionModel(16, 5, 5, 9, [511 << 9], C)
        enum = CEnumerator(model.distances, 9, max_nodes=1000000)
        counts, h_ranks = Counter(), Counter()
        complete = True
        try:
            for xs in enum.leaves():
                nf = normal_form(xs, model.Gamma)
                if any(nf['Gamma'][a] >> b & 1 for a in [8, 9] for b in range(a + 1, 10)):
                    counts['symplectic_Z_rejection'] += 1
                else:
                    rank_h = len(reduced_basis([v >> 8 for v in nf['offsets']]))
                    h_ranks[rank_h] += 1
                    counts['minimum_word_rejection' if rank_h < 2 else 'requires_mixed_test'] += 1
        except EnumerationLimit:
            complete = False
        report = dict(class_index=i, enumeration_complete=complete, leaves=enum.leaf_count,
                      nodes=enum.nodes, counts=dict(counts), H_ranks=dict(h_ranks),
                      status='COMPLETE' if complete else 'LIMIT',
                      excluded=complete and not counts['requires_mixed_test'])
        reports.append(report)
        print(json.dumps(report), flush=True)
    counts = Counter()
    for report in reports:
        counts.update(report['counts'])
    result = dict(status='BOUNDED_BATCH_FINISHED', results=reports, totals=dict(counts),
                  node_limit_per_class=1000000, minimum_stabilizer_weight_assumption=9,
                  seconds=time.monotonic()-started,
                  input_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Counts sampled leaves rejected by minimum-word assumption; '
                        'limited classes are not exclusions; remaining leaves need mixed constraints')
    (BASE / 'minimum_word_benchmark.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='results'}))


if __name__ == '__main__':
    main()
