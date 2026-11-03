"""Necessary kernel/rank/hull screen for any code with the prescribed S_T=R."""
def griesmer_max_dimension(length, distance):
    if distance < 1 or length < 0:
        raise ValueError('Require positive distance and nonnegative length')
    total = 0
    for q in range(1, length + 2):
        total += (distance + (1 << (q - 1)) - 1) // (1 << (q - 1))
        if total > length:
            return q - 1
    raise AssertionError('Unreachable Griesmer sum')


def suffix_stabilizer_bound(length, distance):
    if length < 0 or distance < 1:
        raise ValueError('Invalid suffix stabilizer parameters')
    if distance > length:
        return 0
    # Concatenating each quaternary symbol with (x,z,x+z) doubles its weight.
    # N is a suffix stabilizer, hence an isotropic subspace of dimension <=f.
    bound = min(length, griesmer_max_dimension(3 * length, 2 * distance))
    denominator = 4 * distance - 3 * length
    if denominator > 0:
        plotkin_size = (4 * distance) // denominator
        bound = min(bound, plotkin_size.bit_length() - 1)
    return bound


def kernel_screen(n, k, d, m, r, hull=None, *, minimum_stabilizer_weight=None,
                  localized_low_weights=False):
    if not (0 < k < n and 1 <= m <= n and 1 <= r <= min(m, n - k) and d >= 1):
        raise ValueError('Invalid parameters')
    a, b, cdim = 2 * m - 2 * r, n - k - r, n + k - 2 * m + r
    if minimum_stabilizer_weight is not None and not 1 <= minimum_stabilizer_weight <= m:
        raise ValueError('A full-support local word of weight m bounds the stabilizer minimum')
    suffix_distance = max(d if localized_low_weights else 1, minimum_stabilizer_weight or 1)
    suffix_bound = suffix_stabilizer_bound(n - m, suffix_distance)
    if not 0 <= cdim <= 2 * (n - m):
        return dict(status='EXCLUDED', reason='suffix_dimension', allowed_hulls=[])
    mechanism = None
    if r == 1 and d > m // 2:
        delta = 2 * (d - m // 2)
        bound = griesmer_max_dimension(m, delta)
        pairs = [(q, t) for q in range(bound + 1) for t in range(0, q + 1, 2)]
        mechanism = dict(type='binary_kernel_Griesmer', binary_distance=delta,
                         kernel_dimension_bound=bound)
    elif r == 2 and m % 2 == 0 and d > m - (m + 3) // 4:
        pairs = [(0, 0)]
        mechanism = dict(type='rank_two_cost_radius', kernel_dimension_bound=0)
    elif r == 2 and m % 4 == 0 and d == 3 * (m // 4):
        bound = m.bit_length() - 1
        if (m // 4) % 2:
            bound = min(bound, 2)
            pairs = [(0, 0), (1, 0), (2, 2)]
        else:
            pairs = [(q, 0) for q in range(bound + 1)]
        mechanism = dict(type='balanced_rank_two_kernel', kernel_dimension_bound=bound,
                         odd_symbol_count=bool((m // 4) % 2))
    else:
        return dict(status='UNRESOLVED', reason='no_kernel_bound_in_this_regime')
    allowed = []
    for q, t in pairs:
        rank = a - 2 * q + t
        h = b - rank
        suffix_dimension = b - a + q
        if (q - t // 2 <= k and 0 <= rank <= b and 0 <= h <= min(b, cdim, n - m)
                and 0 <= suffix_dimension <= suffix_bound):
            allowed.append(dict(kernel_dimension=q, kernel_form_rank=t,
                                gram_rank=rank, hull=h,
                                suffix_stabilizer_dimension=suffix_dimension))
    hulls = sorted({p['hull'] for p in allowed})
    excluded = not allowed or (hull is not None and hull not in hulls)
    return dict(status='EXCLUDED' if excluded else 'UNRESOLVED',
                local_rank=r, local_support=m, allowed_hulls=hulls,
                possibilities=allowed, **mechanism,
                suffix_stabilizer_distance=suffix_distance,
                suffix_stabilizer_dimension_bound=suffix_bound,
                minimum_stabilizer_weight_assumption=minimum_stabilizer_weight,
                localized_low_weights_assumption=localized_low_weights,
                scope='Necessary extension conditions with R=S_T; surviving possibilities '
                      'need not be realizable; rank-two R must be full X^m,Z^m')
