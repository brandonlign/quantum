# Reproducing the results

Run every command from the repository root with assertions enabled (never
`python3 -O`). The code checks need only the Python standard library. The
nonexistence proofs use Python 3.12 and the pinned packages:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r verification/requirements.txt
shasum -a 256 -c SHA256SUMS
```

Each command below prints `PASS` (or `REPLAY OK`) and exits with status 0 on
success. Any failed check raises an error. Recorded outputs are in
`results/release_checks/` and `results/recorded_checks/`.

## [[17,5,5]]: an explicit pure code

```sh
python3 codes/17-5-5/verification/verify_short.py
python3 codes/17-5-5/verification/check_paper.py
python3 codes/17-5-5/verification/verify.py
```

`verify_short.py` checks rank 12 and isotropy, finds 1,276 distinct syndromes
for the errors of weight at most two (so the distance is at least five, with
no degeneracy), and checks a weight-five logical operator. `check_paper.py`
reads the generator table, the witness, the upper-bound certificate and the
F4-nonlinearity remark from `papers/sixteen/main.tex` and checks them.
`verify.py` enumerates all 4,194,304 normalizer vectors. Further checks are
listed in [`codes/17-5-5/README.md`](codes/17-5-5/README.md).

## [[18,3,6]]: three explicit pure codes

```sh
python3 codes/18-3-6/verification/verify_all.py
python3 codes/18-3-6/verification/check_paper.py
python3 codes/18-3-6/published-example-audit/audit_cel_published.py
```

`verify_all.py` enumerates all 32,768 stabilizer and 2,097,152 normalizer
vectors of each code and checks the MacWilliams identities, the logical
witnesses, the gluing construction of code A and the transformations giving B
and C. Expected output: `codes/18-3-6/verification/EXPECTED_VERIFICATION.txt`.
`check_paper.py` reads the generator table, enumerator, syndrome counts,
witness and upper-bound certificate from `papers/seventeen/paper/main.tex` and
checks them. The third command checks that the 1999 example of Cohen, Encheva and Litsyn
cannot contain the self-dual subcode its construction needs. See
[`codes/18-3-6/README.md`](codes/18-3-6/README.md) for the C++ error checks.

## [[16,5,5]]: nonexistence

```sh
PYTHON="$PWD/.venv/bin/python" bash verification/release16/replay.sh
```

This takes about 25 minutes on an Apple M3. It checks the Farkas
certificates and the parity identity, runs positive controls on a known
[[16,4,5]] code, runs the complete searches for m = 6 and 7 by two methods, and
spot-checks m = 8 and 9. The options add the rest of the proof:

| Option | Adds | Time |
|---|---|---|
| `--census` | completeness of `census/cells.json` (needs pynauty) | minutes |
| `--full` | the complete m = 8 search and the SAT runs for m = 7, 8 | 17 core-hours |
| `--all` | `--census`, `--full` and both complete m = 9 searches | 45 core-hours |

Every search fails the script unless every class is excluded, and a search
stopped by a node limit counts as not excluded. Logs of the complete runs are
in `verification/release16/search/logs/`.

## [[17,3,6]]: nonexistence

### Finite checks

```sh
.venv/bin/python verification/replay_certificates.py
.venv/bin/python verification/check_constructions.py
.venv/bin/python verification/check_kernel_argument.py
.venv/bin/python verification/check_local_classification.py
.venv/bin/python verification/check_full_plane_orbits.py
.venv/bin/python source/check_record_store.py
```

These take a few seconds. `replay_certificates.py` checks the exact
enumerator certificates (purity and the light-word bound) using only the
standard library. Run the local classification before the orbit check, which
uses its output. Reports go to `audit/`.

### Saved searches

The search records are in the evidence archive on Zenodo,
[doi:10.5281/zenodo.23904905](https://doi.org/10.5281/zenodo.23904905).
Extract both archives in a separate directory and run from its root:

```sh
tar -xzf quantum-review-code.tar.gz
zstd -dc quantum-review-evidence.tar.zst | tar -xf -
cd quantum-review
python3.12 -m venv .venv
.venv/bin/python -m pip install -r verification/requirements.txt
.venv/bin/python -B check_release.py
.venv/bin/python -B source/verify.py
.venv/bin/python -B source/verify_all_record_hashes.py
.venv/bin/python -B source/verify.py --deep
```

`verify.py` re-aggregates every branch of the search from the saved records
and replays the records of eight residual classes (about 70 s).
`verify_all_record_hashes.py` checks the size, CRC and SHA-256 of every
stored record (about 3 minutes). `verify.py --deep` replays every residual
record (about 2.2 hours). Recorded results:
`results/release_checks/saved-replay.json`, `record-integrity.json` and
`deep-replay.json`.

### Fresh searches

The programs that produced the saved records are under
`source/experimenting/qcode-open-entries/general-criterion/unified/`. They
take many core-hours. Fresh runs need their own output paths.

## Papers

```sh
tectonic papers/sixteen/main.tex
bash papers/seventeen/paper/build.sh
```

Built with Tectonic 0.17.0. Read a gzipped report with
`gzip -dc path/to/report.json.gz`.
