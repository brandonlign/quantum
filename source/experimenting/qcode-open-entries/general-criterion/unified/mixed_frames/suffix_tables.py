"""All 56 suffix classes, including hulls 0, 2 and 6; physical replay.

python3 -B -m mixed_frames.suffix_tables --verify-physical
python3 -B -m mixed_frames.suffix_tables --replay
"""
import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path

from model import dual, gram, reduced_basis, span, symp, weight
from .frames import FRAMES, mixed_frame, xor_table

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'suffix9_classes.json'
TABLES = HERE / 'suffix_tables.json'
CERTIFICATE = HERE / 'capacity_certificate.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(C, f=9):
    C = reduced_basis(C)
    P = dual(C, f)
    def syndrome(w):
        return sum(symp(h, w, f) << j for j, h in enumerate(P))
    steps = [(syndrome(w), w) for i in range(f)
             for w in (1 << i, 1 << (f+i), (1 << i) | (1 << (f+i)))]
    delta, reps = [-1] * (1 << len(P)), [0] * (1 << len(P))
    delta[0] = 0
    queue = deque([0])
    while queue:
        s = queue.popleft()
        for t, w in steps:
            u = s ^ t
            if delta[u] < 0:
                delta[u], reps[u] = delta[s]+1, reps[s] ^ w
                queue.append(u)
    assert min(delta) == 0
    Gamma = gram(P, f)
    return dict(C=C, P=P, Gamma=Gamma,
                hull=len(P)-len(reduced_basis(Gamma)), distances=delta,
                representatives=reps, C_distance=min((weight(c, f) for c in span(C) if c), default=None))


def verify_physical(table, f=9):
    """Enumerate every physical word and independently assign its syndrome.

    Also enumerate each representative+C fiber (no BFS). Both minima must
    match. The physical-word oracle checks coverage and the stored lifts.
    """
    P, C = table['P'], table['C']
    elementary = [sum(symp(h, 1 << j, f) << i for i, h in enumerate(P))
                  for j in range(2*f)]
    # Independent elementary commutation spelling: X_j pairs with h's Z_j,
    # Z_j with h's X_j. No symp() in this second calculation.
    explicit = [sum(((h >> (f+j)) & 1) << i for i, h in enumerate(P)) for j in range(f)]
    explicit += [sum(((h >> j) & 1) << i for i, h in enumerate(P)) for j in range(f)]
    assert elementary == explicit
    syndromes = xor_table(elementary)
    physical = [f+1] * (1 << len(P))
    counts = [0] * len(physical)
    mask = (1 << f)-1
    for w, s in enumerate(syndromes):
        physical[s] = min(physical[s], ((w & mask) | (w >> f)).bit_count())
        counts[s] += 1
    assert physical == table['distances']
    assert set(counts) == {1 << len(C)}
    codewords = span(C)
    assert all(syndromes[c] == 0 for c in codewords)
    direct_C_distance = min((sum(bool(((c >> i) | (c >> (f+i))) & 1) for i in range(f))
                             for c in codewords if c), default=None)
    assert direct_C_distance == table['C_distance']
    coset_checks = 0
    for s, rep in enumerate(table['representatives']):
        assert syndromes[rep] == s
        physical_min = f+1
        for c in codewords:
            w = rep ^ c
            physical_min = min(physical_min, ((w & mask) | (w >> f)).bit_count())
            coset_checks += 1
        assert physical_min == table['distances'][s]
    return dict(physical_words_checked=len(syndromes), coset_words_checked=coset_checks)


def capacity_records(tables):
    records = []
    for name in FRAMES:
        frame = mixed_frame(name)
        # Every nonzero Z class costs <6, so X is injective. Distinct cost-one
        # classes need distinct distance-at-least-five syndromes.
        assert max(frame.cost[0]) < 6
        required = sum(c == 1 for c in frame.cost[0])
        for t in tables:
            labels = [s for s, d in enumerate(t['distances']) if d >= 5]
            records.append(dict(frame=name, class_index=t['class_index'],
                                required=required, available=len(labels),
                                distance_five_labels=labels,
                                verdict='PREFILTER_REJECT' if len(labels) < required else 'CANDIDATE'))
    return records


def generate(physical=True):
    source = json.loads(SOURCE.read_text())
    assert len(source['classes']) == 56
    assert sorted(c['class_index'] for c in source['classes']) == list(range(56))
    tables = []
    for c in source['classes']:
        t = build(c['representative']['generators'])
        assert len(t['C']) == 6 and len(t['P']) == 12 and t['C_distance'] >= 6
        assert t['hull'] == c['representative']['hull']
        if physical:
            t['physical_replay'] = verify_physical(t)
        tables.append(dict(class_index=c['class_index'], **t))
    assert Counter(t['hull'] for t in tables) == {0: 15, 2: 33, 6: 8}
    output = dict(status='PHYSICALLY_VERIFIED' if physical else 'UNVERIFIED',
                  source_sha256=sha(SOURCE), label_convention='s_j=<P_j,w>', classes=tables)
    TABLES.write_text(json.dumps(output, indent=1)+'\n')
    certificate = dict(type='DISTANCE_FIVE_CAPACITY', d=6,
                       source_sha256=sha(SOURCE), tables_sha256=sha(TABLES),
                       proof='All nonzero pure-Z costs <6 force injectivity; distinct cost-one classes require distinct syndromes with delta>=5.',
                       records=capacity_records(tables),
                       scope='Necessary prefilter conditional on the supplied 56-class census; candidates are unresolved.')
    CERTIFICATE.write_text(json.dumps(certificate, indent=2)+'\n')
    return output


def replay():
    saved = json.loads(TABLES.read_text())
    cert = json.loads(CERTIFICATE.read_text())
    assert saved['source_sha256'] == cert['source_sha256'] == sha(SOURCE)
    assert cert['type'] == 'DISTANCE_FIVE_CAPACITY' and cert['d'] == 6
    assert cert['tables_sha256'] == sha(TABLES)
    source = json.loads(SOURCE.read_text())['classes']
    assert len(saved['classes']) == len(source) == 56
    for original, t in zip(source, saved['classes']):
        assert t['class_index'] == original['class_index']
        fresh = build(original['representative']['generators'])
        assert all(t[k] == v for k, v in fresh.items())
        assert verify_physical(t) == t['physical_replay']
    for name in FRAMES:
        mixed_frame(name).verify_physical()
    assert cert['records'] == capacity_records(saved['classes'])
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-physical', action='store_true', help='default; always enabled')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    saved = replay() if args.replay else generate()
    records = json.loads(CERTIFICATE.read_text())['records']
    failures = {name: [r['class_index'] for r in records if r['frame'] == name
                       and r['verdict'] == 'PREFILTER_REJECT'] for name in FRAMES}
    summary = dict(status='PASS', replay=args.replay, classes=len(saved['classes']),
                   hull_counts=dict(Counter(t['hull'] for t in saved['classes'])),
                   prefilter_rejects=failures,
                   candidate_pairs=sum(r['verdict'] == 'CANDIDATE' for r in records),
                   physical_words_checked=sum(t['physical_replay']['physical_words_checked'] for t in saved['classes']),
                   coset_words_checked=sum(t['physical_replay']['coset_words_checked'] for t in saved['classes']))
    summary['sha256'] = {p.name: sha(p) for p in [SOURCE, TABLES, CERTIFICATE, Path(__file__),
                                                HERE/'frames.py', HERE.parent/'model.py']}
    (HERE/'suffix_verification.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
