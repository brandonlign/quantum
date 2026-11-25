# Verification

- `release16/`: the [[16,5,5]] proof: Farkas certificates, parity identity,
  lengthening census and extension searches. `release16/replay.sh` runs them.
- `certificates/16/`, `certificates/17/`, `certificates/17-frame-counts/`:
  exact enumerator certificates, checked by `replay_certificates.py`.
- `check_constructions.py`, `check_kernel_argument.py`,
  `check_local_classification.py`, `check_full_plane_orbits.py`: finite checks
  for the [[17,3,6]] proof. Run the local classification before the orbit check.
- `data/`: inputs for those checks.

Commands and expected output are in [REPRODUCING.md](../REPRODUCING.md).
