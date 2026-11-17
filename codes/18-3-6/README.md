# Three explicit pure [[18,3,6]] stabilizer codes

Brandon Li, November 2026.

`manuscript.pdf` is the paper (source `manuscript.tex`, `references.bib`).

## Results

Codes A, B and C each encode three qubits into eighteen with distance exactly
six, matching the known upper bound. A has stabilizer minimum weight 6
(a6 = 4); B and C have stabilizer minimum weight 7 (a7 = 16 and 12). The three
stabilizer enumerators show they are pairwise inequivalent under products of
single-qubit unitaries and qubit permutations. Code A is built by gluing an
8-qubit and a 10-qubit code along a linear quotient map.

The [[18,3,6]] example of Cohen, Encheva and Litsyn (IEEE Trans. Inf. Theory
45, 1999) is invalid: its printed twelve-row extension matrix has maximum
totally isotropic dimension 7, so it cannot contain the self-dual
dimension-9 subcode the construction needs.

## Run the checks

Python 3.10 or newer, standard library only, without `-O`:

```sh
python3 verification/verify_all.py
```

This enumerates all 32,768 stabilizer and 2,097,152 normalizer vectors of each
code, checks the 19 MacWilliams identities, the logical witnesses, the gluing
inequalities and reconstruction, and the transformations relating B and C to A.
Expected output: `verification/EXPECTED_VERIFICATION.txt`. Each of
`verification/frozen_A`, `frozen_B`, `frozen_C` also has a self-contained
`verify_standalone.py` and a C++ check of all errors of weight at most five:

```sh
c++ -O3 -std=c++17 verification/frozen_A/check_errors.cpp -o check_A && ./check_A
```

## Folders

- `verification/`: generators, enumerators and verifiers for A, B and C.
- `reconstruction/`: the gluing construction and search replay for code A, and
  the derivation of B and C (`reconstruction/README.md`).
- `published-example-audit/`: the check of the 1999 printed matrix
  (`python3 published-example-audit/audit_cel_published.py`).
