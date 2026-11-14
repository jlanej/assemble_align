#!/usr/bin/env bash
# Fetch public files listed in a manifest, in order, with resume and size checks.
#
#   fetch_manifest.sh MANIFEST.tsv DATA_ROOT [CONNECTIONS]
#
# MANIFEST.tsv: tab-separated "url <TAB> path relative to DATA_ROOT"; '#' lines are comments.
# Large files are fetched as CONNECTIONS (default 6) byte ranges in parallel: a single long-lived
# connection to these buckets was seen to slow to a fraction of the link. Each range resumes where
# it stopped and is retried when it stalls. A file is renamed into place only when its byte count
# matches the server's Content-Length, so a finished file is always complete. Re-running skips
# finished files. Nothing is written inside the git repository.
set -uo pipefail
manifest=${1:?manifest.tsv}
root=${2:?data root}
nconn=${3:-6}
mkdir -p "$root"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

remote_size() {
  curl -sIL --retry 5 --retry-delay 3 "$1" | tr -d '\r' |
    awk 'tolower($1)=="content-length:"{s=$2} END{print s+0}'
}
local_size() { if [ -e "$1" ]; then stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; else echo 0; fi; }

fetch_range() {  # url, file, first byte, last byte
  local url=$1 file=$2 first=$3 last=$4 want=$(( $4 - $3 + 1 )) have attempt
  for attempt in $(seq 1 40); do
    have=$(local_size "$file")
    [ "$have" -ge "$want" ] && return 0
    # give up on a connection slower than 200 kB/s for 20 s and open a new one
    curl -sfL --speed-limit 200000 --speed-time 20 -r "$(( first + have ))-$last" "$url" >> "$file" || sleep 2
  done
  return 1
}

fail=0
while IFS=$'\t' read -r url rel; do
  case "$url" in ''|\#*) continue ;; esac
  dest="$root/$rel"
  mkdir -p "$(dirname "$dest")"
  if [ -s "$dest" ]; then log "have  $rel"; continue; fi
  want=$(remote_size "$url")
  if [ "$want" -le 0 ]; then log "ERROR no Content-Length for $url"; fail=1; continue; fi
  log "fetch $rel ($(awk -v b="$want" 'BEGIN{printf "%.2f GB", b/1e9}'))"
  n=$nconn; [ "$want" -lt 268435456 ] && n=1
  seg=$(( (want + n - 1) / n ))
  rm -f "$dest.part"
  for i in $(seq 0 $(( n - 1 ))); do
    first=$(( i * seg )); last=$(( first + seg - 1 )); [ "$last" -ge "$want" ] && last=$(( want - 1 ))
    fetch_range "$url" "$dest.seg$i" "$first" "$last" &
  done
  wait
  got=0; for i in $(seq 0 $(( n - 1 ))); do got=$(( got + $(local_size "$dest.seg$i") )); done
  if [ "$got" -eq "$want" ]; then
    : > "$dest.part"
    for i in $(seq 0 $(( n - 1 ))); do cat "$dest.seg$i" >> "$dest.part" && rm -f "$dest.seg$i"; done
    mv "$dest.part" "$dest"; log "done  $rel"
  else
    log "ERROR incomplete $rel ($got of $want bytes); segments kept for the next run"; fail=1
  fi
done < "$manifest"
log "manifest finished (fail=$fail)"
exit $fail
