#!/usr/bin/env bash
# Compare the downloaded inputs with their published checksums (inputs.md5).
#   check_inputs.sh DATA [PATTERN]      PATTERN: only paths that match it (grep -E)
set -uo pipefail
D=${1:?data root}; pat=${2:-.}
here=$(cd "$(dirname "$0")" && pwd)
sum() { if command -v md5sum > /dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
bad=0
while read -r want rel; do
  case "$want" in ''|\#*) continue ;; esac
  printf '%s\n' "$rel" | grep -Eq "$pat" || continue
  if [ ! -s "$D/$rel" ]; then echo "missing   $rel"; bad=1; continue; fi
  got=$(sum "$D/$rel")
  if [ "$got" = "$want" ]; then echo "md5 ok    $rel"; else echo "MISMATCH  $rel: $got, published $want"; bad=1; fi
done < "$here/inputs.md5"
exit $bad
