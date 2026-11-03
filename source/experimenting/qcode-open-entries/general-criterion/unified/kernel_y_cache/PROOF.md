# Cache the immutable local frame

The local R1 frame, basis and physical cost table depend only on m. Neither
the Y solver nor Frame.model_columns mutates them. This copy of the validated
kernel_y_parity solver changes exactly one import: standard(m) is wrapped in
an eight-entry LRU cache. Every search constraint and every solution check is
unchanged. It never caches a search result or a time/node LIMIT. Source
generation asserts the single replacement, and the source digest is pinned.

Repeated leaves for one suffix model therefore avoid reconstructing 2^(2m-2)
physical quotient costs. Complete outputs and node counts must be identical;
time-limited traversal can progress farther and remains unresolved on LIMIT.
