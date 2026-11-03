"""Exact binary algebra for the localized stabilizer extension framework.

Physical words use x | (z << length). Matrices and coordinate vectors use
integer bit masks. This implementation does not import the exploratory filters.
"""
from collections import Counter, deque
from dataclasses import dataclass


def reduced_basis(words):
    piv = {}
    for word in words:
        v = word
        for j in sorted(piv, reverse=True):
            if v >> j & 1:
                v ^= piv[j]
        if v:
            j = v.bit_length() - 1
            for q in list(piv):
                if piv[q] >> j & 1:
                    piv[q] ^= v
            piv[j] = v
    return [piv[j] for j in sorted(piv, reverse=True)]


def span(words):
    result = [0]
    for g in reduced_basis(words):
        result += [v ^ g for v in result]
    return result


def nullspace(rows, width):
    rows = reduced_basis(rows)
    piv = {v.bit_length() - 1: v for v in rows}
    result = []
    for j in range(width):
        if j not in piv:
            v = 1 << j
            for q, row in piv.items():
                if row >> j & 1:
                    v |= 1 << q
            result.append(v)
    return result


def symp(u, v, length):
    mask = (1 << length) - 1
    return (((u & mask) & (v >> length)).bit_count()
            + ((v & mask) & (u >> length)).bit_count()) & 1


def weight(v, length):
    return ((v & ((1 << length) - 1)) | (v >> length)).bit_count()


def dual(words, length):
    mask = (1 << length) - 1
    return nullspace([(v >> length) | ((v & mask) << length)
                      for v in words], 2 * length)


def coordinates(word, words):
    piv = {}
    for i, g in enumerate(words):
        v, c = g, 1 << i
        for j in sorted(piv, reverse=True):
            if v >> j & 1:
                v ^= piv[j][0]
                c ^= piv[j][1]
        if v:
            piv[v.bit_length() - 1] = (v, c)
    v, c = word, 0
    for j in sorted(piv, reverse=True):
        if v >> j & 1:
            v ^= piv[j][0]
            c ^= piv[j][1]
    if v:
        raise ValueError("Word is outside the specified span")
    return c


def lift_coordinate(c, words):
    out = 0
    for i, g in enumerate(words):
        if c >> i & 1:
            out ^= g
    return out


def embed(prefix, suffix, m, f):
    x = (prefix & ((1 << m) - 1)) | ((suffix & ((1 << f) - 1)) << m)
    z = (prefix >> m) | ((suffix >> f) << m)
    return x | (z << (m + f))


def split(word, m, f):
    n = m + f
    x, z = word & ((1 << n) - 1), word >> n
    p = (x & ((1 << m) - 1)) | ((z & ((1 << m) - 1)) << m)
    w = (x >> m) | ((z >> m) << f)
    return p, w


def quotient_basis(R, m):
    """Physical lifts of a symplectic basis of R^perp/R, Z halves first."""
    R = reduced_basis(R)
    assert all(symp(u, v, m) == 0 for u in R for v in R)
    available, seen = [], list(R)
    for v in dual(R, m):
        if len(reduced_basis(seen + [v])) > len(seen):
            available.append(v)
            seen.append(v)
    zs, xs = [], []
    while available:
        u = available.pop(0)
        i = next(i for i, v in enumerate(available) if symp(u, v, m))
        v = available.pop(i)
        zs.append(u)
        xs.append(v)
        available = [w ^ (u if symp(w, v, m) else 0)
                       ^ (v if symp(w, u, m) else 0) for w in available]
    out = zs + xs
    half = len(zs)
    assert len(out) == 2 * m - 2 * len(R)
    assert all(symp(u, v, m) == int((i < half and j == i + half)
                                   or (j < half and i == j + half))
               for i, u in enumerate(out) for j, v in enumerate(out))
    return out


def gram(words, length):
    return [sum(symp(u, v, length) << j for j, v in enumerate(words))
            for u in words]


def pullback_gram(columns, b):
    half = len(columns) // 2
    out = [0] * b
    for x, y in zip(columns[:half], columns[half:]):
        for j in range(b):
            if x >> j & 1:
                out[j] ^= y
            if y >> j & 1:
                out[j] ^= x
    return out


@dataclass
class ExtensionModel:
    n: int
    k: int
    d: int
    m: int
    R: list
    C: list

    def __post_init__(self):
        assert 0 <= self.k <= self.n and 0 <= self.m <= self.n and self.d >= 1
        self.f = self.n - self.m
        assert all(0 <= v < (1 << (2 * self.m)) for v in self.R)
        assert all(0 <= v < (1 << (2 * self.f)) for v in self.C)
        self.R = reduced_basis(self.R)
        self.C = reduced_basis(self.C)
        self.r = len(self.R)
        assert len(self.C) == self.n + self.k - 2 * self.m + self.r
        self.E = quotient_basis(self.R, self.m)
        self.a = len(self.E)
        self.P = dual(self.C, self.f)
        self.b = len(self.P)
        assert self.b == self.n - self.k - self.r
        self.Gamma = gram(self.P, self.f)
        self.hull = self.b - len(reduced_basis(self.Gamma))
        self.cost = [min(weight(lift_coordinate(e, self.E) ^ r, self.m)
                         for r in span(self.R)) for e in range(1 << self.a)]
        self.C_distance = min((weight(v, self.f) for v in span(self.C) if v),
                              default=None)
        self.distances, self.representatives = self._coset_distances()

    def syndrome(self, w):
        return sum(symp(h, w, self.f) << j for j, h in enumerate(self.P))

    def _coset_distances(self):
        N = 1 << self.b
        distances, representatives = [-1] * N, [0] * N
        distances[0] = 0
        steps = []
        for i in range(self.f):
            for word in [1 << i, 1 << (self.f + i),
                         (1 << i) | (1 << (self.f + i))]:
                steps.append((self.syndrome(word), word))
        queue = deque([0])
        while queue:
            s = queue.popleft()
            for t, w in steps:
                u = s ^ t
                if distances[u] == -1:
                    distances[u] = distances[s] + 1
                    representatives[u] = representatives[s] ^ w
                    queue.append(u)
        assert all(d >= 0 for d in distances)
        assert all(self.syndrome(w) == s for s, w in enumerate(representatives))
        return distances, representatives

    def map_value(self, e, columns):
        return lift_coordinate(e, columns)

    def criterion(self, columns):
        assert len(columns) == self.a
        assert all(0 <= c < (1 << self.b) for c in columns)
        suffix = self.C_distance is None or self.C_distance >= self.d
        compatibility = pullback_gram(columns, self.b) == self.Gamma
        weights = all(self.cost[e] + self.distances[self.map_value(e, columns)] >= self.d
                      for e in range(1, 1 << self.a))
        return suffix and compatibility and weights

    def construct(self, columns):
        D = [embed(r, 0, self.m, self.f) for r in self.R]
        D += [embed(0, c, self.m, self.f) for c in self.C]
        D += [embed(e, self.representatives[s], self.m, self.f)
              for e, s in zip(self.E, columns)]
        D = reduced_basis(D)
        assert len(D) == self.n + self.k
        return dual(D, self.n), D

    def direct_property(self, columns):
        S, D = self.construct(columns)
        isotropic = all(symp(u, v, self.n) == 0 for u in S for v in S)
        Rglobal = set(embed(r, 0, self.m, self.f) for r in span(self.R))
        outside_R = all(weight(v, self.n) >= self.d
                        for v in span(D) if v not in Rglobal)
        ST = {v for v in span(S) if split(v, self.m, self.f)[1] == 0}
        assert ST == Rglobal
        return isotropic and outside_R

    def kernel_rank_identity(self, columns):
        equations = [sum(((c >> j) & 1) << i for i, c in enumerate(columns))
                     for j in range(self.b)]
        K = nullspace(equations, self.a)
        q = len(K)
        Kphysical = [lift_coordinate(e, self.E) for e in K]
        t = len(reduced_basis(gram(Kphysical, self.m)))
        rank_G = len(reduced_basis(pullback_gram(columns, self.b)))
        assert rank_G == self.a - 2 * q + t
        return q, t, rank_G

    def small_prefix_obstruction(self):
        """Necessary hull/capacity conditions when max prefix cost < d."""
        assert max(self.cost) < self.d
        required_hull = self.b - self.a
        if self.hull != required_hull:
            return {"excluded": True, "reason": "hull", "hull": self.hull,
                    "required_hull": required_hull}
        Q = span(self.Gamma)  # im Gamma = N^perp/C when K=0 and N=rad(P)
        assert len(Q) == 1 << self.a
        counts = {}
        for j in range(self.d):
            actual = sum(self.distances[s] <= j for s in Q if s)
            capacity = sum(c >= self.d - j for c in self.cost[1:])
            counts[j] = {"actual": actual, "capacity": capacity}
            if actual > capacity:
                return {"excluded": True, "reason": "capacity", "threshold": j,
                        "counts": counts}
        return {"excluded": False, "reason": "passes_necessary_tests", "counts": counts}

    def prefix_profile(self):
        return dict(sorted(Counter(self.cost).items()))


def recover_extension(S, n, k, d, m):
    """Recover R, C and the exact map from a stabilizer matrix, without a
    centralizer enumeration. The selected sites must already be first.
    """
    S = reduced_basis(S)
    assert len(S) == n - k
    assert all(symp(u, v, n) == 0 for u in S for v in S)
    f = n - m
    projected = [split(v, m, f)[1] for v in S]
    equations = [sum(((w >> j) & 1) << i for i, w in enumerate(projected))
                 for j in range(2 * f)]
    R = [split(lift_coordinate(c, S), m, f)[0]
         for c in nullspace(equations, len(S))]
    C = dual(reduced_basis(projected), f)
    model = ExtensionModel(n, k, d, m, R, C)
    graph = {}
    for v in dual(S, n):
        p, w = split(v, m, f)
        e = coordinates(p, model.R + model.E) >> model.r
        s = model.syndrome(w)
        for j in sorted(graph, reverse=True):
            if e >> j & 1:
                e ^= graph[j][0]
                s ^= graph[j][1]
        if e:
            graph[e.bit_length() - 1] = (e, s)
        else:
            assert s == 0
    columns = []
    for i in range(model.a):
        e, s = 1 << i, 0
        for j in sorted(graph, reverse=True):
            if e >> j & 1:
                e ^= graph[j][0]
                s ^= graph[j][1]
        assert e == 0
        columns.append(s)
    assert set(span(model.construct(columns)[0])) == set(span(S))
    return model, columns
