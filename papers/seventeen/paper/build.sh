#!/usr/bin/env bash
# Build main.pdf with Tectonic (https://tectonic-typesetting.github.io).
set -euo pipefail
cd "$(dirname -- "$0")"
# Fix the PDF metadata timestamp so the same source and compiler reproduce it.
export SOURCE_DATE_EPOCH=1794787200
"${TECTONIC:-tectonic}" main.tex
