# Preserved source and saved evidence

`verify.py` runs the mathematical evidence checker through `record_store.py`.
The companion archive places the database at `source/evidence/search_records.sqlite`
and supplies the remaining inputs under `source/experimenting/`.

The adapter relocates historical paths and reads original payload bytes.
Programs, inputs and fingerprinted producer binaries retain their provenance;
saved replay does not execute those binaries. Release results are in
`../results/release_checks/`.

See [REPRODUCING.md](../REPRODUCING.md) for commands and scope.
