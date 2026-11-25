# A pure [[17,5,5]] stabilizer code

The code and its proofs are in Section 8 of
[`papers/sixteen/main.pdf`](../../papers/sixteen/main.pdf). The twelve generators
are in `verification/generators.txt` (qubit 0 on the left).

## Checks

Python 3.10 or newer, standard library only, without `-O`. Each runs in a few
seconds from the repository root.

| Command | What it checks |
|---|---|
| `python3 codes/17-5-5/verification/verify_short.py` | rank 12, isotropy, 1,276 distinct radius-two syndromes, a weight-five logical operator: pure distance exactly 5 |
| `python3 codes/17-5-5/verification/check_paper.py` | the generator table, witness, upper-bound certificate and F4-nonlinearity remark as printed in the paper |
| `python3 codes/17-5-5/verification/verify.py` | all 4,194,304 normalizer vectors and all errors of weight at most five; weight enumerators |
| `python3 codes/17-5-5/verification/verify_local_unitary.py` | no local-unitary equivalence to a GF(4)-linear code |
| `python3 codes/17-5-5/verification/construct_a5_orbit.py` | the code is spanned by the orbit of one weight-eight word under a monomial A5 subgroup of order 60 |
| `python3 codes/17-5-5/verification/normal_form/verify_standalone.py` | the A5 presentation, built from its printed cycles alone, gives a pure [[17,5,5]] code |
| `python3 codes/17-5-5/verification/normal_form/verify_equivalence.py` | the A5 presentation is equivalent to the generator table |
| `python3 codes/17-5-5/verification/build_logical_basis.py` | five logical X/Z pairs (`logical_basis.json`) |

Expected output is in `verification/VERIFICATION.json` and `verification/verify.log`.
`documentation/normal_form_verification.txt` describes the A5 presentation.
