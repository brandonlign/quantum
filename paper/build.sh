#!/usr/bin/env bash
set -euo pipefail
paper_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
# Fix the PDF metadata timestamp so the same source and compiler reproduce it.
export SOURCE_DATE_EPOCH=1793491200
tectonic --keep-logs --keep-intermediates "$paper_dir/main.tex"
