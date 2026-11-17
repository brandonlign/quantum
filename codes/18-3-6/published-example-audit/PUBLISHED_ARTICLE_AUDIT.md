# Final published article: independently checked

**Finding: the final 1999 IEEE article retains the same claim and the same obstructed matrix. The displayed example was not repaired between the inspected preprint and the final journal publication.**

## Primary source

Gérard Cohen, Sylvia Encheva, and Simon Litsyn, “On Binary Constructions of Quantum Codes,” *IEEE Transactions on Information Theory* 45(7), 2495–2498 (November 1999), DOI [10.1109/18.796389](https://doi.org/10.1109/18.796389).

The final IEEE PDF is publicly hosted by the authors' institution:
https://biblio.telecom-paristech.fr/cgi-bin/download.cgi?id=9892

The inspected PDF has four pages, 128,437 bytes, and SHA-256:
`9e9f47adf47e625563ca7f2c3f279d6ceb47a1cdf55a05c22ce36acc14f67395`.

The journal page header, IEEE production metadata, article title, and author list identify this as the final publication, not the arXiv preprint. The article begins partway through page 2495, following the end of another article. The relevant Theorem 1, Table I, Example 1, and matrix are all on printed page 2496 (second PDF page, zero-based index 1).

## What was directly inspected

The final PDF page was rendered and independently read. Theorem 1 still requires C-perpendicular contained in C. Example 1 still explicitly claims a [[18,3,6]] code from a binary [18,9,6] code C extended to a [18,12,4] code C-prime. Table I repeats these parameters. All twelve binary matrix rows are identical to the previously audited preprint rows; this equality was checked visually and by exact text extraction from the identified final PDF.

## Recomputed result and proof

The first nine displayed rows span a [18,9,6] binary code but have Euclidean Gram rank 8. Rows 1 and 2 overlap in positions 3, 5, and 6, so they have inner product 1 and are not self-orthogonal.

More strongly, the full printed twelve-row matrix spans a [18,12,4] binary code. Its restricted Euclidean form is alternating, has rank 10, and has a 2-dimensional radical. A totally isotropic subspace has at most dimension 5 after passage to the nondegenerate 10-dimensional quotient, plus at most dimension 2 in the radical. Thus every totally isotropic subspace of the displayed C-prime has dimension at most 7.

A length-18 dimension-9 binary C with C-perpendicular contained in C must be self-dual, hence totally isotropic of dimension 9. No such subcode can lie inside this displayed C-prime. This obstruction does not depend on assigning the first nine rows to C, on the choice of an enlargement map, or on the ordering of coordinates.

Consequently, **the data printed in the final published Example 1 do not satisfy the hypothesis of the theorem invoked to justify that example**.

## Exact executable audit

Run the standalone calculation with no external dependencies:

    python3 audit_cel_published.py

To additionally verify the final PDF hash and extract its twelve printed rows, with Poppler installed:

    python3 audit_cel_published.py /path/to/Cohen_Encheva_Litsyn_1999_institution9892.pdf

The second form was run successfully during review. `CEL_PUBLISHED_MATRIX_AUDIT.json` and `audit_cel_published.log` include the complete matrix, both complete binary weight enumerators, Gram data, radical certificates, source identification, and the exact comparison with the preprint audit.

## Consequences and limits

The previous uncertainty about whether the final journal version fixes the displayed example is resolved: it does not. References or manuscript text saying the final article was unavailable or uninspected should be updated.

This finding concerns the specific displayed construction. The positive certificates for the present three codes, and their product-unitary/permutation inequivalence, are independent of it.
