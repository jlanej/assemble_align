#!/usr/bin/env bash
# Build the report PDF from markdown: pandoc (md -> standalone HTML) then weasyprint (HTML -> PDF).
# Usage: ./build_pdf.sh report.md report.pdf   (run inside an env with pandoc + weasyprint)
set -euo pipefail
in=${1:-assemble_align_review.md}; out=${2:-${in%.md}.pdf}
here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/aa_pdf.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
pandoc "$in" -f markdown+pipe_tables+grid_tables+implicit_figures+raw_html+smart -t html5 --standalone \
  --toc --toc-depth=2 --number-sections --css "$here/style.css" --embed-resources --resource-path="$(dirname "$in"):$here" \
  --lua-filter "$here/locallinks.lua" --metadata lang=en -o "$tmp/report.html"
weasyprint --base-url "$(dirname "$in")" "$tmp/report.html" "$out" 2> >(grep -v -E "^WARNING: (Ignored|Expected)" >&2 || true)
ls -l "$out" | awk '{printf "PDF bytes: %s\n", $5}'
