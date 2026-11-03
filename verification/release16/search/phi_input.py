"""Write phisearch input for a given dmin list and m."""
from full import prefix
pop = lambda v: bin(v).count('1')
def make_input(dmin, m, d=5):
    basis, cost = prefix(m)
    nb = len(basis)
    low = []
    for c in range(1, 1 << nb):
        v = 0
        for j in range(nb):
            if c >> j & 1: v ^= basis[j]
        t = d - cost(v)
        if t > 0: low.append((c, t))
    s = f'{m} {len(dmin)}\n' + ' '.join(map(str, dmin)) + f'\n{len(low)}\n' + '\n'.join(f'{c} {t}' for c, t in low) + '\n'
    return s
