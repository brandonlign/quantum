# Quantum stabilizer codes

Brandon Li. This repository settles four open entries in Grassl's table of
binary stabilizer codes ([codetables.de](https://codetables.de/QECC/QECC.php?q=4)).

| Parameters | Result | Table before |
|---|---|---|
| [[16,5,5]] | does not exist, including degenerate codes; d_max(16,5) = 4 | 4–5 |
| [[17,5,5]] | explicit pure code; d_max(17,5) = 5 | 4–5 |
| [[17,3,6]] | does not exist, including degenerate codes; d_max(17,3) = 5 | 5–6 |
| [[18,3,6]] | three explicit pure codes; d_max(18,3) = 6 | 5–6 |

## Papers

- [`papers/sixteen/main.pdf`](papers/sixteen/main.pdf): *The shortest binary
  stabilizer code with five logical qubits and distance five*. Proves that no
  [[16,5,5]] code exists and gives the [[17,5,5]] code.
- [`papers/seventeen/paper/main.pdf`](papers/seventeen/paper/main.pdf): proves
  that no [[17,3,6]] code exists.
- [`codes/18-3-6/manuscript.pdf`](codes/18-3-6/manuscript.pdf): the three
  [[18,3,6]] codes.

## Checking each result

Run from the repository root with assertions enabled (no `python3 -O`).

| Result | Command | Time | Needs |
|---|---|---|---|
| [[17,5,5]] | `python3 codes/17-5-5/verification/verify_short.py` | 1 s | standard library |
| [[18,3,6]] | `python3 codes/18-3-6/verification/verify_all.py` | 3 s | standard library |
| [[16,5,5]] | `bash verification/release16/replay.sh` | 25 min | numpy, sympy, python-sat, C compiler |
| [[16,5,5]], whole proof | `bash verification/release16/replay.sh --all` | 45 core-hours | as above, plus pynauty |
| [[17,3,6]], finite checks | `python3 verification/replay_certificates.py` and the four `verification/check_*.py` | 2 s | numpy, pynauty |
| [[17,3,6]], saved searches | `python source/verify.py` | 70 s | evidence archive from Zenodo |
| [[17,3,6]], every saved record | `python source/verify.py --deep` | 2.2 h | evidence archive from Zenodo |

[REPRODUCING.md](REPRODUCING.md) gives the setup, every command and the
expected output. Recorded outputs are in `results/`.

## Layout

| Folder | Contents |
|---|---|
| `papers/sixteen/` | the [[16,5,5]]/[[17,5,5]] paper (LaTeX source and PDF) |
| `papers/seventeen/` | the [[17,3,6]] paper |
| `codes/17-5-5/` | generators and checks for the [[17,5,5]] code |
| `codes/18-3-6/` | generators, checks and construction of the [[18,3,6]] codes |
| `verification/release16/` | the [[16,5,5]] certificates, census and searches (`replay.sh`) |
| `verification/` | certificates and finite checks for the [[17,3,6]] proof |
| `source/` | the [[17,3,6]] search programs, saved evidence and its checker |
| `results/` | recorded outputs of the checks; large reports are gzipped |

The full evidence archive for the [[17,3,6]] searches is on Zenodo,
[doi:10.5281/zenodo.23904905](https://doi.org/10.5281/zenodo.23904905).

Papers and data: CC BY 4.0. Code: MIT. See [LICENSES.md](LICENSES.md).
