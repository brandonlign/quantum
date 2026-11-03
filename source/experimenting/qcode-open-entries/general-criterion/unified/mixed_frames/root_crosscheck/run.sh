#!/bin/bash
# Root cross-check: stage-1 Python reference (mu=6, NO mu_out prune) on 668 classes excluded by stage 2.
cd "$(dirname "$0")/../.."
for c in "$@"; do
  PYTHONDONTWRITEBYTECODE=1 python3 -B -m mixed_frames.enumerate_z --frame 668 --class-index $c \
    --seconds 3000 --nodes 100000000 --leaves 1000 --output mixed_frames/root_crosscheck/668_$c.json
  python3 -c "import json;d=json.load(open('mixed_frames/root_crosscheck/668_$c.json'));print($c,d['status'],d['nodes'],d['leaf_count'],round(d['seconds'],1))"
done
