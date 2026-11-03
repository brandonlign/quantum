# Run instructions

Use Python 3.12 with assertions enabled (no `-O`). Run from the repository root.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r verification/requirements.txt
shasum -a 256 -c SHA256SUMS
```

## Quick checks

```sh
.venv/bin/python -B verification/replay_certificates.py
.venv/bin/python -B verification/check_constructions.py
.venv/bin/python -B verification/check_kernel_argument.py
.venv/bin/python -B verification/check_local_classification.py
.venv/bin/python -B verification/check_full_plane_orbits.py
.venv/bin/python -B source/check_record_store.py
```

Expect PASS. Run classification before orbits. Reports go to `audit/` and
`source/provenance/`. The first three checks use only the standard library.

## Full evidence

Download the code and evidence archives from
[Zenodo](https://doi.org/10.5281/zenodo.23897639). Extract them together in a
separate directory and run these commands from that release's root:

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

The database belongs at `source/evidence/search_records.sqlite`. These commands
validate saved computations. Deep replay can take several hours.
Original archive inventories and uncompressed reports are included in the release.

## Fresh searches

Requires a C compiler. These runs can take many core-hours.

```sh
PYTHON="$PWD/.venv/bin/python" bash verification/release16/replay.sh --census --full
(cd verification/release16/search && PHIPROG=./phisearch ../../../.venv/bin/python run_c.py 9)
```

The first command builds the engines and checks m=7,8; m=9 runs separately.
Seventeen-qubit producers are under
`source/experimenting/qcode-open-entries/general-criterion/unified/`.
New runs require separate output paths, rebuilt engines and new fingerprints.

## Paper and reports

Build the paper with Tectonic 0.17.0:

```sh
bash paper/build.sh
```

Read a compressed report with `gzip -dc path/to/report.json.gz`.
Preserved source, inputs and fingerprinted engines retain their original bytes.
