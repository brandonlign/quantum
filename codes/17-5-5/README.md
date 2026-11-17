# An optimal pure [[17,5,5]] stabilizer code

Brandon Li, November 2026.

`manuscript/main.pdf` is the paper (LaTeX source `manuscript/main.tex`).
`verification/` contains the explicit code and its checks.

## Results

- Twelve commuting independent generators define five logical qubits.
- Minimum stabilizer weight 8; exact pure distance 5.
- 1,276 distinct radius-two syndromes show that no nonzero normalizer word has
  weight at most four; a displayed weight-five logical word gives equality.
- A rational MacWilliams identity excludes [[17,5,d>=6]] stabilizer codes,
  including degenerate codes, so the code is optimal.
- A three-site odd-overlap witness excludes local-Clifford equivalence to any
  GF(4)-linear code, even allowing permutations; a reduced-projector purity
  certificate extends this to arbitrary local unitaries.
- The code is the span of the orbit of one weight-eight word under a binary
  monomial A5 subgroup of order 60.

## Run the checks

Python 3.10 or newer, standard library only, without `-O`:

```sh
python3 verification/verify_short.py
python3 verification/check_manuscript.py
python3 verification/verify_local_unitary.py
python3 verification/construct_a5_orbit.py
python3 verification/normal_form/verify_standalone.py
python3 verification/normal_form/verify_equivalence.py
python3 verification/verify.py
```

`verify.py` enumerates all 4,194,304 normalizer vectors and all physical errors
through weight five. Expected results are in `verification/VERIFICATION.json`
and `verification/verify.log`. `verification/logical_basis.json` gives five
logical X/Z pairs (`build_logical_basis.py` rebuilds and checks them).
`documentation/normal_form_verification.txt` explains the normal-form
certificates.
