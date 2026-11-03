# Unary consequences of the full kernel equations

This is a separate experimental namespace; frozen production is untouched.

In fixed-Z normal coordinates the low Y columns U_i satisfy
U_j[i]+U_i[j]=G[i,j]. A kernel equation is XOR_{j:beta_j=1} U_j = rhs.
Taking coordinate i and substituting the pair equations gives

    beta dot U_i = rhs_i + (G beta)_i.

For any sector alpha, linearity therefore gives

    beta dot U_alpha = alpha dot (rhs + G beta).

Filter candidate U_alpha by these parity equations BEFORE constructing its
distance domain. Every full Gram-compatible solution with the fixed kernel
satisfies them, so this preserves the complete solution set. No phase gauge
or quotient is used. Existing vector kernel checks, pair checks, all sector
weight checks and full model.criterion assertions remain in place.

The new constraints use only the top-left p by p Gram block and low rhs.
An inconsistent high rhs is rejected by the existing kernel_projection test.
If a kernel has q independent X projections, these q scalar parities reduce
each raw candidate sector from 2^p to 2^(p-q), before weight evaluation.

This improves Y completion only. It does not help a Z tree that has no leaves,
and it does not establish any target exclusion without a new complete census.
