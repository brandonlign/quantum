"""Stage-2 physical graph constraints and frozen enumerator inputs."""
import hashlib
import json
from pathlib import Path

from model import ExtensionModel, lift_coordinate, coordinates, nullspace, weight
from mixed_frames.frames import mixed_frame, xor_table
from mixed_frames.suffix_tables import SOURCE

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARY = HERE / 'enumerate_z'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes():
    paths = [SOURCE, ROOT/'model.py', *HERE.parent.glob('*.py'),
             *HERE.glob('*.py'), HERE/'enumerate_z.cpp', BINARY,
             HERE.parent/'capacity_certificate.json', HERE.parent/'suffix_tables.json',
             ROOT/'weight8_split/frame_count_certificates.json']
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def target(name, index):
    frame = mixed_frame(name)
    c = next(c for c in json.loads(SOURCE.read_text())['classes'] if c['class_index'] == index)
    return ExtensionModel(17, 3, 6, 8, frame.R, c['representative']['generators']), frame


def graph_data(model):
    return [weight(h, model.f) for h in xor_table(model.P)]


def nonsingular_graph_weights(model):
    if model.hull or model.a != model.b:
        return None
    out = [0]*(1 << model.b)
    for a,w in enumerate(graph_data(model)):
        out[lift_coordinate(a,model.Gamma)] = w
    return out


def graph_minimum(model, frame, columns):
    """Exact minimum of S outside R. Coefficient zero labels R itself.

    Prefix adjoint e(a) = (Y^T a, X^T a) in (z,x) order.
    P has independent physical generators even when its Gram form is singular.
    """
    p = frame.p
    weights = graph_data(model)
    return min((frame.costs[sum(((a & c).bit_count() & 1) << ((i+p) % (2*p))
                              for i, c in enumerate(columns))]+weights[a]
                for a in range(1, 1 << model.b)), default=model.n+1)


def graph_z_ok(model, frame, xs, mu):
    if not mu:
        return True
    weights = graph_data(model)
    for a in xor_table(nullspace(xs, model.b))[1:]:
        try:
            z = coordinates(lift_coordinate(a, model.Gamma), xs)
        except ValueError:
            return False
        if frame.costs[z]+weights[a] < mu:
            return False
    return True


def graph_partial_ok(model, frame, xs, mu):
    if not mu:
        return True
    j, Q = len(xs), 1 << frame.p
    caps = [max(frame.costs[z | (x << frame.p)] for x in range(Q)
                if x & ((1 << j)-1) == low for z in range(Q)) for low in range(1 << j)]
    weights = graph_data(model)
    lowcap = min(caps)
    for a in range(1, 1 << model.b):
        if weights[a]+lowcap >= mu:
            continue
        low = sum(((a & c).bit_count() & 1) << i for i, c in enumerate(xs))
        if weights[a]+caps[low] < mu:
            return False
    return True


def input_text(model, frame, mu=0, rank=True, domains=True):
    assert max(frame.cost[0]) < model.d or model.d < 6  # controls may fix X explicitly
    arrays = [[frame.p, model.b, model.d, mu, int(rank), int(domains)],
              frame.costs, model.distances, model.Gamma, graph_data(model)]
    return '\n'.join(' '.join(map(str, a)) for a in arrays)+'\n'


def write_input(path, model, frame, mu=0, rank=True, domains=True):
    Path(path).write_text(input_text(model,frame,mu,rank,domains))


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def prepare():
    dest = HERE/'inputs'
    dest.mkdir(exist_ok=True)
    manifest = []
    records = json.loads((HERE.parent/'capacity_certificate.json').read_text())['records']
    for name in ('668', '778', '688'):
        for r in records:
            if r['frame'] != name or r['verdict'] != 'CANDIDATE':
                continue
            model, frame = target(name, r['class_index'])
            mu = 8 if name in ('668', '778') else 6
            path = dest/f'{name}_{r["class_index"]}.txt'
            write_input(path, model, frame, mu)
            first = [s for s in range(1 << model.b)
                     if frame.costs[1]+model.distances[s] >= model.d]
            manifest.append(dict(frame=name, class_index=r['class_index'], d=6, mu_out=mu,
                                 input=str(path.relative_to(HERE)), input_sha256=sha(path),
                                 first_images=first, R=frame.R, E=frame.E, C=model.C,
                                 P=model.P, Gamma=model.Gamma, hull=model.hull))
    atomic_json(HERE/'manifest.json', dict(classes=manifest, hashes=hashes(),
                label_convention='s_j=<P_j,w>; e=z|(x<<6)'))
    return manifest


if __name__ == '__main__':
    print(json.dumps(dict(prepared=len(prepare()))))
