"""Streaming adapter for the independently checked canonical C enumerator."""
import json
from pathlib import Path
import subprocess
import tempfile
from z_enumerator import EnumerationLimit


class CEnumerator:
    def __init__(self, distances, m, *, seconds=None, max_nodes=None):
        self.distances = distances
        self.m = m
        self.seconds = seconds
        self.max_nodes = max_nodes
        self.nodes = self.leaf_count = 0
        self.process = None
        self.summary = None

    def leaves(self):
        command = [str(Path(__file__).with_name('z_enumerator_c'))]
        if self.seconds is not None:
            command += ['--seconds', str(self.seconds)]
        if self.max_nodes is not None:
            command += ['--max-nodes', str(self.max_nodes)]
        with tempfile.TemporaryFile(mode='w+') as errors:
            self.process = subprocess.Popen(command, stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE, stderr=errors,
                                            text=True, bufsize=1)
            try:
                self.process.stdin.write(f'{self.m} {len(self.distances)} 5\n')
                self.process.stdin.write(' '.join(map(str, self.distances)) + '\n')
                self.process.stdin.close()
                for line in self.process.stdout:
                    fields = line.split()
                    if not fields or fields[0] != 'LEAF' or len(fields) != self.m:
                        raise RuntimeError(f'Invalid C output: {line!r}')
                    self.leaf_count += 1
                    yield tuple(map(int, fields[1:]))
                code = self.process.wait()
                errors.seek(0)
                output = errors.read().strip()
                self.summary = json.loads(output)
                if code != 0 or self.summary['status'] == 'ERROR':
                    raise RuntimeError(f'C enumerator failed: {output}')
                self.nodes = self.summary['nodes']
                if self.summary['leaves'] != self.leaf_count:
                    raise RuntimeError('C/Python emitted-leaf counts disagree')
                if self.summary['status'] == 'LIMIT':
                    raise EnumerationLimit(self.summary['reason'])
                if self.summary['status'] != 'COMPLETE':
                    raise RuntimeError(f'Unrecognized C status: {output}')
            finally:
                if self.process.poll() is None:
                    self.process.terminate()
                    self.process.wait()
                self.process.stdout.close()

