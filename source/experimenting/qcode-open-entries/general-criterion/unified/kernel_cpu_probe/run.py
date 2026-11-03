"""Bounded same-source M3/PGO compiler experiment; production is untouched."""
import hashlib
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha


def call(binary, case, *, profile_path=None, reps=None):
    model, data, context, chosen = case['input']
    command = [str(binary), '--gamma', '--avoid-radical', '--first-orbits', '--kernel-data']
    if case['backend'] == 'joint':
        command += ['--mitm4', '65536']
    command += ['--seconds', '12']
    env = dict(os.environ)
    if profile_path is not None:
        env['LLVM_PROFILE_FILE'] = str(profile_path)
    with tempfile.TemporaryDirectory() as temporary:
        group = Path(temporary) / 'group.txt'
        group.write_text(f'{len(context["group"])} 7\n' +
                         '\n'.join(' '.join(map(str, a)) for a in context['group']) + '\n')
        command += ['--prefix-group', str(group)]
        started = time.monotonic()
        p = subprocess.run(command, input=base.payload(model, data, chosen if reps is None else reps, context),
                           env=env, text=True, capture_output=True, check=True, timeout=35)
        wall = time.monotonic() - started
    stats = json.loads(p.stderr)
    return dict(wall=wall, stats=stats, stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest())


def main():
    assert __debug__
    plan = base.validate(require_controls=True)
    rows = {r['cls']: r for r in plan['records']}
    sources = {'compact': ROOT / 'kernel_z_compact/z_enum_compact.c',
               'joint': ROOT / 'kernel_z_joint/z_enum_joint.c'}
    originals = {'compact': base.Z_BINARY, 'joint': base.Z_JOINT}
    frozen = {str(p.relative_to(ROOT)): sha(p) for p in list(sources.values()) + list(originals.values())}
    profdata = subprocess.check_output(['xcrun', '--find', 'llvm-profdata'], text=True).strip()
    compiler = subprocess.check_output(['cc', '--version'], text=True)
    commands = []
    for backend, source in sources.items():
        for mode in ('native', 'instrumented'):
            binary = HERE / f'{mode}_{backend}'
            command = ['cc', '-O3', '-std=c11', '-mcpu=native', '-fomit-frame-pointer', '-flto']
            if mode == 'instrumented':
                command += ['-fprofile-instr-generate']
            command += [str(source), '-o', str(binary)]
            subprocess.run(command, check=True)
            commands.append(command)
    selected = [(0, p) for p in (0, 1, 3, 8, 12)] + [(5, p) for p in (4, 5, 13)]
    cases = []
    for cls, profile in selected:
        row = rows[cls]
        s = next(s for s in expected_slices(plan, row, profiles=[profile]) if s['oi'] == 0)
        model, data, H, order, contexts, io = base.slice_setup(tuple(row['C']), tuple(s['branch']['N_orbits'][0]['N']),
                            tuple(s['prefix']['H']), s['prefix']['p'], s['q'], s['t'])
        context = contexts[s['ci']]
        ss = base.suffix_slice(model, s['branch']['N_orbits'][0]['N'], data, context)
        reps = ss['first_representatives'][:64]
        cases.append(dict(cls=cls, profile=profile, backend='joint' if profile in (0, 8) else 'compact',
                          input=(model, data, context, reps)))
    for i, case in enumerate(cases):
        path = HERE / f'train_{case["backend"]}_{i}.profraw'
        result = call(HERE / f'instrumented_{case["backend"]}', case, profile_path=path,
                      reps=case['input'][3][:2])
        print(json.dumps(dict(stage='training', cls=case['cls'], profile=case['profile'], status=result['stats']['status'])), flush=True)
    for backend, source in sources.items():
        raw = sorted(HERE.glob(f'train_{backend}_*.profraw'))
        merged = HERE / f'{backend}.profdata'
        subprocess.run([profdata, 'merge', '-output', str(merged), *map(str, raw)], check=True)
        command = ['cc', '-O3', '-std=c11', '-mcpu=native', '-fomit-frame-pointer', '-flto',
                   '-fprofile-instr-use=' + str(merged), str(source), '-o', str(HERE / f'pgo_{backend}')]
        subprocess.run(command, check=True)
        commands.append(command)
    records = []
    for round_number in (0, 1):
        for case in cases:
            modes = ('baseline', 'native', 'pgo') if round_number == 0 else ('pgo', 'native', 'baseline')
            results = {}
            for mode in modes:
                binary = originals[case['backend']] if mode == 'baseline' else HERE / f'{mode}_{case["backend"]}'
                results[mode] = call(binary, case)
            if all(x['stats']['status'] == 'COMPLETE' for x in results.values()):
                reference = results['baseline']
                for mode, result in results.items():
                    assert result['stdout_sha256'] == reference['stdout_sha256'], (case['cls'], case['profile'], mode, 'stream')
                    assert {k: v for k, v in result['stats'].items() if k != 'seconds'} == {
                        k: v for k, v in reference['stats'].items() if k != 'seconds'}, (case['cls'], case['profile'], mode, 'counters')
            record = dict(round=round_number, cls=case['cls'], profile=case['profile'], backend=case['backend'], results=results)
            records.append(record)
            print(json.dumps(dict(stage='comparison', **{k: v for k, v in record.items() if k != 'results'},
                                  walls={k: round(v['wall'], 4) for k, v in results.items()},
                                  status={k: v['stats']['status'] for k, v in results.items()})), flush=True)
            atomic_json(HERE / 'benchmark.json', dict(status='BOUNDED_COMPILER_PROBE_IN_PROGRESS', records=records,
                                                    source_hashes=frozen, scope='No production change or target exclusion'))
    assert frozen == {str(p.relative_to(ROOT)): sha(p) for p in list(sources.values()) + list(originals.values())}
    complete = [r for r in records if all(v['stats']['status'] == 'COMPLETE' for v in r['results'].values())]
    summary = {}
    for mode in ('native', 'pgo'):
        ratios = [r['results'][mode]['wall'] / r['results']['baseline']['wall'] for r in complete]
        summary[mode] = dict(median_ratio=statistics.median(ratios),
                             total_ratio=sum(r['results'][mode]['wall'] for r in complete) /
                                         sum(r['results']['baseline']['wall'] for r in complete))
    atomic_json(HERE / 'benchmark.json', dict(status='BOUNDED_COMPILER_PROBE_COMPLETE', records=records,
                complete_comparisons=len(complete), summary=summary, commands=commands, compiler=compiler,
                source_hashes=frozen, hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__)] +
                    [HERE / f'{mode}_{backend}' for mode in ('native', 'pgo') for backend in sources]},
                scope='Same source compiler experiment with paired complete stream/counter checks. No production deployment or target exclusion.'))
    print(json.dumps(dict(stage='summary', complete_comparisons=len(complete), summary=summary)), flush=True)


if __name__ == '__main__':
    main()
