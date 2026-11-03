# Prefix data are independent of the suffix Gram matrix

In fixed_kernel.common.kernel_contexts the input Gamma is used only to check
its dimension, check its radical dimension, and set rad=nullspace(Gamma,b).
All remaining calculations use H,m: K and its radical, anchor site order,
normalized H, the GL(g,2) images, tags, projected prefix Gram matrices, physical
automorphisms and closed inner groups. None of those uses Gamma or rad.

Compute that prefix-only data once using a canonical alternating Gamma with
the required dimension and radical dimension. For an actual input, reproduce
the original dimension/radical checks and substitute precisely its original
nullspace(Gamma,b), in precisely the same order, into every context. Everything
else is unchanged. Enumeration order and every returned value match the
original function. Cached group lists are read-only in all callers.

The outer identification-orbit function uses only H,m and the ordered tuple
of radical_images from each context. Its exact result can therefore be cached
by those inputs too. Return a copy so a caller cannot mutate its saved result.
No search outcome, syndrome distance, first-label orbit or LIMIT is cached.

This changes setup cost only. It cannot exclude an extension or alter a search
tree. Frozen fixed_kernel sources remain untouched.
