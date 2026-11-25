# Three pure [[18,3,6]] stabilizer codes

Code A and its proofs are in Section 10 of
[`papers/seventeen/paper/main.pdf`](../../papers/seventeen/paper/main.pdf), with
codes B and C in its Remark 10.4. Each code encodes three qubits into eighteen
with distance exactly six, which is optimal. A has stabilizer minimum weight 6
(4 words of weight 6); B and C have stabilizer minimum weight 7 (16 and 12
words of weight 7). Their stabilizer weight enumerators differ, so no two of
them are related by qubit permutations and single-qubit unitaries.

## Checks

Python 3.10 or newer, standard library only, without `-O`, from the repository
root:

| Command | What it checks | Time |
|---|---|---|
| `python3 codes/18-3-6/verification/verify_all.py` | all 32,768 stabilizer and 2,097,152 normalizer vectors of each code, the MacWilliams identities, logical witnesses, the gluing construction of A, and the transformations giving B and C | 3 s |
| `python3 codes/18-3-6/verification/check_paper.py` | the generator table, enumerator, syndrome counts, witness and upper-bound certificate as printed in the paper | 1 s |
| `python3 codes/18-3-6/published-example-audit/audit_cel_published.py` | the printed 1999 matrix of Cohen, Encheva and Litsyn has no self-dual [18,9] subcode | 1 s |

Expected output of `verify_all.py`: `verification/EXPECTED_VERIFICATION.txt`.
Each of `verification/frozen_A`, `frozen_B` and `frozen_C` also has a
self-contained `verify_standalone.py` and a C++ check of all errors of weight at
most five:

```sh
c++ -O3 -std=c++17 codes/18-3-6/verification/frozen_A/check_errors.cpp -o check_A && ./check_A
```

## Folders

- `verification/`: generators, enumerators and checks for A, B and C.
- `reconstruction/`: the gluing construction and search for code A, and the
  derivation of B and C (`reconstruction/README.md`).
- `published-example-audit/`: the check of the 1999 example
  ([`PUBLISHED_ARTICLE_AUDIT.md`](published-example-audit/PUBLISHED_ARTICLE_AUDIT.md)).
