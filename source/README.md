# [[17,3,6]] search programs and saved evidence

- `experimenting/qcode-open-entries/general-criterion/unified/`: the programs
  that ran the [[17,3,6]] searches, their inputs and their aggregated results.
  Files keep their original bytes, because the checker verifies their
  fingerprints.
- `quantum-extension/verification/verify_evidence.py`: the evidence checker.
- `verify.py`: runs the checker. `record_store.py` serves the saved search
  records from the evidence database and maps the original absolute paths to
  this folder.
- `verify_all_record_hashes.py`: checks every stored record against its hash.
- `benchmarks/`: timing comparisons, with results in `../results/benchmarks/`.

The search records themselves are in the Zenodo evidence archive
([doi:10.5281/zenodo.23904905](https://doi.org/10.5281/zenodo.23904905)), which
places the database at `source/evidence/search_records.sqlite`. Commands are in
[REPRODUCING.md](../REPRODUCING.md).
