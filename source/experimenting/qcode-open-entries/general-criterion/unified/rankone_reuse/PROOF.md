# Exact reuse of the q1 profile-1 enumeration for profile 2

This version changes enumeration cost, not mathematical coverage. It imports the
frozen residual runner and strict chunk validator without changing their files.
The only reused tasks are q=1,t=0, identification 0, profile 2. Their source is
profile 1 of the same suffix class, N orbit, and first-label chunk. Each source
must pass the original strict chunk validator as EXCLUDED_SLICE_CHUNK, use the
compact binary, have known source fingerprints, and contain its entire COMPLETE
Z stream. No LIMIT, FOUND, incomplete Y evidence, or reuse chain is admissible.

Before reuse the implementation reconstructs both task contexts from the frozen
plan, including suffix symmetries and first-label chunks. It compares the actual
serialized stdin and prefix-group bytes, binary hash and flags. A hash alone
never establishes input equality. The native executable is deterministic for
these unlimited invocations: its only inputs are argv, stdin, and the prefix
group file; its clock changes only timing without an enabled time limit. Thus
identical inputs have identical ordered leaves and counters except timing.
The independent bounded paired runs in opus_parallel/m9_acceleration/zreuse/
support this argument; they are supplementary evidence, not an extra axiom.

Reuse returns precisely the source Z statistics and ordered X leaves to the
frozen run_task. That function still solves each leaf using the target's own H,
persists raw FOUND data, checks physical witnesses, and assembles target evidence.
It is never permissible to copy a source Y rejection. A process-local temporary
replacement of z_search is restored in a finally block; each pool worker runs
one task at a time. The source record is hashed before and after use. Its path,
hash, input/group hashes, binary hash, argv and reconstructed stdout hash are
recorded. z.seconds is explicitly inherited; target seconds is actual elapsed
task time including source verification. The immutable plan and baseline source
fingerprints are cached per worker after process admission; the actual source
record and executable bytes are still hashed for each reuse. The executable
digest must also equal the digest in the source producer's fingerprints.
This version makes no claim of a new physical positive control at
distance six, whose existence is the unresolved research question.

New fingerprints include the engine source, proof, verification source and PASS
report. They differ from frozen legacy fingerprints. Pool provenance separately
pins pool, records, aggregate and their integration verification. Existing
legacy records are never rewritten. The new aggregator reconstructs all original
expected tasks, invokes the original strict validator for each record, and also
reconstructs and revalidates every reuse source. It checks ordered leaf equality
and source digests again. Unknown producers, conflicting completed records,
missing tasks, LIMIT, malformed records, and saved FOUND data keep coverage open
or erroneous. In particular an orphaned candidate cannot become an exclusion.

The coverage theorem, classifications, physical symmetries, and the other q1/q2
and q3 branches retain their existing proof obligations. Reuse alone excludes
no additional mathematical case.
