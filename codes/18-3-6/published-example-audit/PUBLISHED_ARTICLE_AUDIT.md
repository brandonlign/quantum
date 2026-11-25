# The 1999 [[18,3,6]] example of Cohen, Encheva and Litsyn

Gérard Cohen, Sylvia Encheva and Simon Litsyn, "On Binary Constructions of
Quantum Codes," *IEEE Transactions on Information Theory* 45(7), 2495–2498
(1999), DOI [10.1109/18.796389](https://doi.org/10.1109/18.796389). The final
PDF is hosted by the authors' institution at
https://biblio.telecom-paristech.fr/cgi-bin/download.cgi?id=9892 (four pages,
128,437 bytes, SHA-256
`9e9f47adf47e625563ca7f2c3f279d6ceb47a1cdf55a05c22ce36acc14f67395`).

Example 1 on page 2496 claims an [[18,3,6]] code from a binary [18,9,6] code C
inside an [18,12,4] code C′, using Theorem 1, which requires C⊥ ⊆ C. The twelve
printed rows of C′ are the same in the arXiv preprint (quant-ph/9812065) and
the journal version.

## Why the printed data do not give the code

The twelve printed rows span an [18,12,4] binary code. The ordinary dot product
restricted to it is alternating, with rank 10 and a two-dimensional radical, so
every self-orthogonal subcode has dimension at most 2 + 10/2 = 7. A
length-18, dimension-9 code C with C⊥ ⊆ C is self-dual, hence self-orthogonal
of dimension 9, so no such C lies inside the printed C′. This holds for every
choice of C, enlargement map and coordinate order. (The first nine rows alone
have Gram rank 8: rows 1 and 2 overlap in positions 3, 5 and 6.)

The positive results for codes A, B and C do not depend on this example.

## Commands

    python3 audit_cel_published.py
    python3 audit_cel_published.py /path/to/final.pdf

The first recomputes everything from the printed rows using only the standard
library. The second also checks the PDF hash and extracts the rows (needs
Poppler's `pdftotext`). Output: `CEL_PUBLISHED_MATRIX_AUDIT.json` and
`audit_cel_published.log`.
