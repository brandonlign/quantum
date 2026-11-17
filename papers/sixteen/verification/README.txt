Verification for "No binary [[16,5,5]] stabilizer code exists"

cells.json
  Representatives of the shortened codes: keys "10,10", "9,8", "8,6", "7,4"
  for m = 6, 7, 8, 9 (37, 2298, 1579, 43 classes).
  SHA-256 181b27547b89c25c0ff783dc8ac9cb926e4d8b2d78e9631a8854780d4d114755

The certificates, census and extension searches are in
verification/release16/. From that folder:

  ./replay.sh            certificates, parity identity, positive controls,
                         complete m=6,7 searches, m=8 and m=9 spot checks
                         (about 25 minutes on an Apple M3); ends "REPLAY OK"
  ./replay.sh --census   also re-checks completeness of cells.json (needs pynauty)
  ./replay.sh --full     also reruns the complete m=8 search and the m=7,8 SAT runs
                         (about 17 core-hours)
  cd search && PHIPROG=./phisearch python3 run_c.py 9
                         complete m=9 search (about 14 core-hours)

Logs of the original complete runs are in search/logs/ (see README.md there).
