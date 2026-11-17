#!/usr/bin/env bash
set -euo pipefail
paper_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
# Fix the PDF metadata timestamp so the same source and compiler reproduce it.
export SOURCE_DATE_EPOCH=1794787200
cd "$paper_dir"
/opt/homebrew/bin/tectonic --only-cached --keep-logs --keep-intermediates main.tex
