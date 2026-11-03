"""Reuse an identical COMPLETE native stream; solve Y separately for each H."""
import copy
import hashlib
import json
import time
from functools import lru_cache
from pathlib import Path

from rankone_residual import common as base
from rankone_residual.aggregate import verify_chunk
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ENGINE_FILES = ['rankone_reuse/common.py', 'rankone_reuse/PROOF.md',
                'rankone_reuse/verify_engine.py']


def source_hashes():
    return dict(base.fingerprints(), **{name: sha(ROOT / name) for name in ENGINE_FILES})


@lru_cache(maxsize=1)
def baseline_hashes():
    return base.fingerprints()


@lru_cache(maxsize=1)
def saved_plan():
    # Validated at process admission and immutable during a production sweep.
    return json.loads(base.PLAN.read_text())


def fingerprints():
    report_path = HERE / 'engine_verification.json'
    report = json.loads(report_path.read_text())
    assert report['status'] == 'PASS' and report['hashes'] == source_hashes()
    return dict(source_hashes(), **{'rankone_reuse/engine_verification.json': sha(report_path)})


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def native_input(model, data, context, reps):
    """The exact frozen compact invocation, with only the temp path normalized."""
    return dict(binary=str(base.Z_BINARY.relative_to(ROOT)), binary_sha256=sha(base.Z_BINARY),
                argv=['--gamma', '--avoid-radical', '--first-orbits', '--kernel-data',
                      '--prefix-group', '<group bytes>'],
                stdin=base.payload(model, data, reps, context).encode(),
                group=(f'{len(context["group"])} {model.m-1}\n' +
                       '\n'.join(' '.join(map(str, a)) for a in context['group']) + '\n').encode())


def setup(row, s, k):
    prefix, branch = s['prefix'], s['branch']
    N = branch['N_orbits'][s['oi']]['N']
    model, data, H, order, contexts, io = base.slice_setup(
        tuple(row['C']), tuple(N), tuple(prefix['H']), prefix['p'], s['q'], s['t'])
    context = contexts[s['ci']]
    ss = base.suffix_slice(model, N, data, context)
    mine = ss['first_representatives'][k * base.CHUNK:(k + 1) * base.CHUNK]
    assert 0 <= k < s['chunks'] and base.digest(ss) == s['orbit_digest']
    return model, data, context, mine


def verify_source(path, row, s, k, accepted_hashes):
    path = Path(path).resolve()
    path.relative_to(ROOT)
    assert (s['q'], s['t'], s['prefix']['profile_index'], s['ci']) == (1, 0, 2, 0)
    plan = saved_plan()
    sources = [x for x in expected_slices(plan, row, branches={(1, 0)}, profiles=[1])
               if (x['oi'], x['ci']) == (s['oi'], 0)]
    assert len(sources) == 1
    source_slice = sources[0]
    assert path.name == source_slice['names'][k]
    before = sha(path)
    r = json.loads(path.read_text())
    assert r['hashes'] in accepted_hashes, 'source producer fingerprints'
    assert 'z_reuse' not in r, 'reuse chains are prohibited'
    assert r['backend'] == 'compact'
    if 'reuse_pool_hashes' in r:
        from rankone_reuse.records import pool_hashes
        assert r['reuse_pool_hashes'] == pool_hashes(), 'source pool fingerprints'
    assert verify_chunk(path, row, source_slice, k, r['hashes']) == 'EXCLUDED_SLICE_CHUNK'
    left = native_input(*setup(row, source_slice, k))
    right = native_input(*setup(row, s, k))
    assert left == right, 'native inputs must match byte for byte'
    assert left['binary_sha256'] == r['hashes']['kernel_z_compact/z_enum_compact'], 'source binary changed'
    leaves = [tuple(e['X']) for e in r['leaves']]
    stdout = ''.join('LEAF ' + ' '.join(map(str, xs)) + '\n' for xs in leaves).encode()
    assert sha(path) == before, 'source changed during verification'
    metadata = dict(source_path=str(path.relative_to(ROOT)), source_sha256=before,
                    binary_sha256=left['binary_sha256'], argv=left['argv'],
                    input_sha256=digest_bytes(left['stdin']), group_sha256=digest_bytes(left['group']),
                    stdout_sha256=digest_bytes(stdout),
                    timing='z.seconds inherited from source; seconds measures this target task')
    return r, leaves, left, metadata


def run_task(row, s, k, hashes, source_path=None, candidate_dir=None, backend='compact'):
    start = time.monotonic()
    args = (row, s['branch'], s['oi'], s['prefix'], s['ci'], k, None, hashes)
    if source_path is None:
        r, _ = base.run_task(*args, candidate_dir=candidate_dir, backend=backend)
        return r
    assert backend == 'compact'
    source, leaves, native, metadata = verify_source(source_path, row, s, k,
                                                    [baseline_hashes(), hashes])
    original = base.z_search

    def replay(model, data, c, reps, seconds, *, group=True, fixed=None, max_nodes=None, backend='compact'):
        assert seconds is None and group and fixed is None and max_nodes is None and backend == 'compact'
        assert native_input(model, data, c, reps) == native
        return copy.deepcopy(source['z']), list(leaves)

    try:
        base.z_search = replay
        r, _ = base.run_task(*args, candidate_dir=candidate_dir, backend=backend)
    finally:
        base.z_search = original
    assert sha(Path(source_path)) == metadata['source_sha256'], 'source changed during Y solve'
    r['z_reuse'] = metadata
    r['seconds'] = time.monotonic() - start
    # The frozen runner persists a raw candidate before physical verification.
    # Add origin metadata to that already-safe record without losing its columns.
    if 'candidate_file' in r:
        path = ROOT / r['candidate_file']
        raw = json.loads(path.read_text())
        raw['record'] = r
        atomic_json(path, raw)
    return r


def validate_reuse(path, row, s, k, accepted_hashes):
    r = json.loads(Path(path).read_text())
    assert r['backend'] == 'compact'
    metadata = r['z_reuse']
    source, leaves, native, expected = verify_source(ROOT / metadata['source_path'], row, s, k, accepted_hashes)
    assert metadata == expected, 'reuse origin metadata'
    assert r['z'] == source['z'], 'reused Z evidence and timing'
    xs = [tuple(e['X']) for e in r['leaves']]
    if r['status'] == 'EXCLUDED_SLICE_CHUNK':
        assert xs == leaves, 'target must solve every source leaf'
    else:
        assert xs == leaves[:len(xs)], 'target leaves must be an ordered source prefix'
    return expected
