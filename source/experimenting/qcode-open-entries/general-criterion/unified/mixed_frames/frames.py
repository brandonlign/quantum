"""Adapted physical quotient bases. Coordinate e = z | (x << p).

Run from the repository root: python3 -B -m mixed_frames.frames.
The oracle enumerates physical centralizer words, rather than R translates.
"""
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path

from model import (coordinates, dual, gram, lift_coordinate, nullspace,
                   reduced_basis, span, symp, weight)

FRAMES = {'668': [0xff00, 0x300f], '688': [0xff00, 0x003f],
          '778': [0xff00, 0x403f]}
HISTOGRAMS = {'668': [1, 4, 64, 360, 1175, 1556, 840, 96],
              '688': [1, 2, 74, 354, 1285, 1660, 720],
              '778': [1, 2, 70, 378, 1355, 1930, 360]}


def xor_table(columns):
    out = [0]
    for c in columns:
        out += [v ^ c for v in out]
    return out


@dataclass
class Frame:
    m: int
    R: list
    E: list
    costs: list
    name: str = 'control'

    @property
    def p(self):
        return len(self.E) // 2

    @property
    def cost(self):
        q = 1 << self.p
        return [self.costs[x*q:(x+1)*q] for x in range(q)]

    def model_columns(self, model, columns):
        """Change domain coordinates, leaving suffix syndrome labels unchanged."""
        assert model.m == self.m and model.R == self.R
        return [lift_coordinate(coordinates(v, self.R + self.E) >> len(self.R),
                                columns) for v in model.E]

    def adapted_columns(self, model, columns):
        assert model.m == self.m and model.R == self.R
        return [lift_coordinate(coordinates(v, model.R + model.E) >> model.r,
                                columns) for v in self.E]

    def verify_physical(self):
        """All 4^m physical words; quotient coordinates from pairings only."""
        q = 1 << self.p
        minima = [self.m + 1] * (q*q)
        counts = [0] * (q*q)
        for w in range(1 << (2*self.m)):
            if any(symp(w, r, self.m) for r in self.R):
                continue
            z = sum(symp(w, f, self.m) << i
                    for i, f in enumerate(self.E[self.p:]))
            x = sum(symp(w, a, self.m) << i
                    for i, a in enumerate(self.E[:self.p]))
            e = z | (x << self.p)
            # Independent physical weight spelling.
            wt = sum(bool(((w >> i) | (w >> (self.m+i))) & 1)
                     for i in range(self.m))
            minima[e] = min(minima[e], wt)
            counts[e] += 1
        assert minima == self.costs
        assert set(counts) == {1 << len(self.R)}
        return sum(counts)

    def record(self):
        return dict(name=self.name, m=self.m, R=self.R, E=self.E,
                    convention='e=z|(x<<p); E=Z+F; cost[x][z]',
                    costs=self.costs, histogram=dict(Counter(self.costs)),
                    pure_z_histogram=dict(Counter(self.cost[0])))


def adapted_frame(m, R, name='control'):
    R = reduced_basis(R)
    assert all(symp(u, v, m) == 0 for u in R for v in R)
    p = m - len(R)
    Z, seen = [], list(R)
    for z in nullspace([r & ((1 << m)-1) for r in R], m):
        w = z << m
        if len(reduced_basis(seen + [w])) > len(seen):
            Z.append(w)
            seen.append(w)
    if len(Z) != p:
        raise ValueError('The pure-Z quotient is not a full Lagrangian')
    # Fixed labels, ordered once by their physical singleton costs. In the
    # target frames this places a cost-one vector first without sorting images.
    translates = span(R)
    Z.sort(key=lambda w: (min(weight(w ^ r, m) for r in translates), w))
    available = span(dual(R, m))
    F = []
    for i in range(p):
        F.append(next(w for w in available
                      if all(symp(z, w, m) == int(i == j)
                             for j, z in enumerate(Z))
                      and all(symp(f, w, m) == 0 for f in F)))
    E = Z + F
    assert gram(E, m) == ([1 << (i+p) for i in range(p)]
                          + [1 << i for i in range(p)])
    costs = [min(weight(w ^ r, m) for r in translates) for w in xor_table(E)]
    return Frame(m, R, E, costs, name)


def mixed_frame(name):
    return adapted_frame(8, FRAMES[name], name)


def main():
    records = []
    for name in FRAMES:
        frame = mixed_frame(name)
        checked = frame.verify_physical()
        h = Counter(frame.costs)
        assert [h[i] for i in range(max(h)+1)] == HISTOGRAMS[name]
        assert max(frame.cost[0]) == 4
        records.append(dict(**frame.record(), physical_words_checked=checked))
    path = Path(__file__).with_name('prefix_tables.json')
    path.write_text(json.dumps(dict(status='PASS', frames=records), indent=2)+'\n')
    print(json.dumps(dict(status='PASS', frames=len(records),
                         centralizer_words_checked=sum(r['physical_words_checked'] for r in records))))


if __name__ == '__main__':
    main()
