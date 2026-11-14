#!/usr/bin/env bash
# Run several fetch manifests one after another (bandwidth is the limit, so never in parallel).
#   fetch_queue.sh DATA_ROOT [--after LOGFILE] MANIFEST...
# With --after, wait until LOGFILE reports that an earlier manifest finished.
set -uo pipefail
root=${1:?data root}; shift
here=$(cd "$(dirname "$0")" && pwd)
if [ "${1:-}" = "--after" ]; then
  until grep -q "manifest finished" "$2" 2>/dev/null; do sleep 20; done
  shift 2
fi
mkdir -p "$root/logs"
for m in "$@"; do
  name=$(basename "$m" .tsv)
  "$here/fetch_manifest.sh" "$m" "$root" > "$root/logs/fetch_$name.log" 2>&1 || true
done
