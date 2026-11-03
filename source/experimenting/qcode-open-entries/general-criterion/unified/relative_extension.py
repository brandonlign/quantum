"""Reference exact logical-distance profile, including nonlocalized impurity."""
from math import inf
from model import lift_coordinate, pullback_gram, weight


def logical_distance(model, columns, *, max_suffix_sites=8):
    if model.f > max_suffix_sites:
        raise ValueError('Reference profile exceeds explicit suffix enumeration limit')
    if len(columns) != model.a or any(not 0 <= c < 1 << model.b for c in columns):
        raise ValueError('Invalid prefix map')
    if pullback_gram(columns, model.b) != model.Gamma:
        return dict(compatible=False, distance=None, profile=None)
    half = model.a // 2
    P = [0] * (1 << model.b)
    labels = {}
    for c in range(1 << model.b):
        if c:
            low = c & -c
            P[c] = P[c ^ low] ^ model.P[low.bit_length()-1]
        adjoint = sum(((c & v).bit_count() & 1) << j for j, v in enumerate(columns))
        tau = (adjoint >> half) | ((adjoint & ((1 << half)-1)) << half)
        labels[P[c]] = tau
    # Store only the two cheapest distinct labels for each syndrome. Excluding
    # one stabilizer-prefix label can then never require more than the runner-up.
    minima = [[] for _ in range(1 << model.b)]
    for w in range(1 << (2 * model.f)):
        s, key, wt = model.syndrome(w), labels.get(w, None), weight(w, model.f)
        entries = minima[s]
        for i, (value, label) in enumerate(entries):
            if label == key:
                entries[i] = (min(value, wt), label)
                break
        else:
            entries.append((wt, key))
        entries.sort(key=lambda item: item[0])
        del entries[2:]
    profile = []
    distance = inf
    for e in range(1 << model.a):
        s = lift_coordinate(e, columns)
        delta = next((value for value, label in minima[s] if label != e), inf)
        profile.append(delta)
        distance = min(distance, model.cost[e] + delta)
    return dict(compatible=True, distance=distance, profile=profile)
