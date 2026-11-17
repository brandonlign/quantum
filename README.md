# Quantum stabilizer codes

Brandon Li. Four open entries in Grassl's table of binary stabilizer codes are settled:

| Parameters | Result | Table before |
|---|---|---|
| [[16,5,5]] | does not exist (including degenerate codes), so d_max(16,5) = 4 | 4–5 |
| [[17,3,6]] | does not exist (including degenerate codes), so d_max(17,3) = 5 | 5–6 |
| [[17,5,5]] | explicit pure code, optimal | 4–5 |
| [[18,3,6]] | three explicit pure codes, optimal | 5–6 |

## Papers

- [`papers/sixteen/main.pdf`](papers/sixteen/main.pdf): No binary [[16,5,5]] stabilizer code exists.
- [`papers/seventeen/paper/main.pdf`](papers/seventeen/paper/main.pdf): No binary [[17,3,6]] stabilizer code exists.
- [`codes/17-5-5/manuscript/main.pdf`](codes/17-5-5/manuscript/main.pdf): an optimal pure [[17,5,5]] code.
- [`codes/18-3-6/manuscript.pdf`](codes/18-3-6/manuscript.pdf): three pure [[18,3,6]] codes.
- [`paper/main.pdf`](paper/main.pdf): the combined [[16,5,5]]/[[17,3,6]] paper ([Zenodo](https://doi.org/10.5281/zenodo.23874369)).

## Folders

- `papers/sixteen/verification/`: the shortened-code classes for the [[16,5,5]] proof; the searches are in `verification/release16/` (`replay.sh`).
- `papers/seventeen/computation/`: the exact rank-one checker and its controls for the [[17,3,6]] proof.
- `codes/17-5-5/`, `codes/18-3-6/`: generators and standalone verifiers.
- `paper/`: the combined paper, PDF and LaTeX source.
- `verification/`: certificates and standalone checks for the nonexistence proofs ([REPRODUCING.md](REPRODUCING.md)).
- `source/`: research programs, inputs and evidence reader.
- `results/`: saved reports and benchmarks; large JSON reports use gzip.

The two constructions check in a few seconds with the Python standard library:

```sh
python3 -B codes/17-5-5/verification/verify_short.py
python3 -B codes/18-3-6/verification/verify_all.py
```

Paper and data: CC BY 4.0. Code: MIT. See [LICENSES.md](LICENSES.md).
