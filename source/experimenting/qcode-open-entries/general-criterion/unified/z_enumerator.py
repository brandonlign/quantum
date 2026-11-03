"""Exact Z-subset enumeration using cached binary syndrome masks.

For R=<Z^m>, d>=3 and m>=3, actual compatible maps can be normalized so
ALL m Z-images are increasing, including the dependent final image. This
canonical mode is complete for actual extensions, but need not preserve the
number of leaves of the exploratory enumerator, which orders only m-1 images.
"""
import time


class EnumerationLimit(RuntimeError):
    pass


class ZEnumerator:
    def __init__(self, distances, m, d=5, *, canonical=True, seconds=None,
                 max_nodes=None):
        self.distances = list(map(int, distances))
        self.N = len(self.distances)
        assert self.N > 0 and self.N & (self.N - 1) == 0
        self.b = self.N.bit_length() - 1
        assert m >= 3 and d >= 3
        self.m, self.d, self.canonical = m, d, canonical
        self.seconds, self.max_nodes = seconds, max_nodes
        self.nodes = self.leaf_count = 0
        self.start = None
        self.all_bits = (1 << self.N) - 1
        self.swap_masks = []
        for j in range(self.b):
            step = 1 << j
            low = sum(((1 << step) - 1) << k for k in range(0, self.N, 2 * step))
            self.swap_masks.append((step, low, self.all_bits ^ low))
        # Bit t of translated[need][s] is 1 exactly when dist[t XOR s]>=need.
        self.translated = {}
        for need in range(1, d):
            good = sum(1 << s for s, dist in enumerate(self.distances) if dist >= need)
            self.translated[need] = [self.xor_permute(good, s) for s in range(self.N)]

    def xor_permute(self, bits, syndrome):
        for j, (step, low, high) in enumerate(self.swap_masks):
            if syndrome >> j & 1:
                bits = ((bits & low) << step) | ((bits & high) >> step)
        return bits

    def leaves(self):
        self.start = time.monotonic()
        if max(self.distances) < self.d - 1:
            return
        if sum(v >= self.d - 1 for v in self.distances) < self.m:
            return

        def rec(xs, sums, sizes, lo):
            self.nodes += 1
            if self.max_nodes is not None and self.nodes > self.max_nodes:
                raise EnumerationLimit('node_limit')
            if self.seconds is not None and self.nodes % 1024 == 0:
                if time.monotonic() - self.start >= self.seconds:
                    raise EnumerationLimit('time_limit')
            allowed = self.all_bits ^ ((1 << (lo + 1)) - 1)
            for s, size in zip(sums, sizes):
                need = self.d - min(size + 1, self.m - size - 1)
                if need > 0:
                    allowed &= self.translated[need][s]
                    if not allowed:
                        return
            while allowed:
                bit = allowed & -allowed
                x = bit.bit_length() - 1
                allowed ^= bit
                if len(xs) == self.m - 2:
                    last = x ^ sums[-1]  # XOR of all previous images plus x.
                    if self.canonical and last <= x:
                        continue
                    self.leaf_count += 1
                    yield xs + [x]
                else:
                    yield from rec(xs + [x], sums + [s ^ x for s in sums],
                                   sizes + [s + 1 for s in sizes], x)

        yield from rec([], [0], [0], -1)
