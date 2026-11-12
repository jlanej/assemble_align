#!/usr/bin/env python3
"""Extract all rows of one query sample from a remote, BGZF-compressed PAF whose rows are grouped
by query sample in sorted order (true for the HPRC release-2 impg all-vs-1 *.merged.paf.gz files).
Binary-searches BGZF block boundaries with HTTP range requests and downloads only the byte range
that holds the sample (~30 MB instead of 6.5 GB).

usage: l5_fetch_paf_range.py URL SAMPLE OUT.paf.gz
"""
import sys, zlib, gzip, urllib.request

def http_size(url):
    req = urllib.request.Request(url, method='HEAD')
    return int(urllib.request.urlopen(req).headers['Content-Length'])

def get(url, a, b):
    req = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (a, b)})
    return urllib.request.urlopen(req).read()

def block_at(url, off, size, n=200000):
    """first complete BGZF block at or after `off`: returns (abs_offset, first_sample)."""
    b = get(url, off, min(size - 1, off + n - 1))
    i = b.find(b'\x1f\x8b\x08\x04')
    while i >= 0:
        try:
            d = zlib.decompressobj(31)
            out = d.decompress(b[i:])
            lines = out.split(b'\n')
            # first complete line (skip a possibly partial first line unless block is at 0)
            cand = lines[0:] if off + i == 0 else lines[1:]
            for l in cand[:-1]:
                if l:
                    return off + i, l.split(b'\t')[0].split(b'#')[0].decode()
            return off + i, None
        except Exception:
            pass
        i = b.find(b'\x1f\x8b\x08\x04', i + 1)
    return None, None

def main(url, sample, out):
    size = http_size(url)
    # lower bound: largest block start whose first sample < sample
    lo, hi = 0, size
    while hi - lo > 200000:
        mid = (lo + hi) // 2
        o, s = block_at(url, mid, size)
        if o is None or s is None or s >= sample:
            hi = mid
        else:
            lo = mid
    start, _ = block_at(url, lo, size)
    # upper bound: smallest block start whose first sample > sample
    lo2, hi2 = start, size
    while hi2 - lo2 > 200000:
        mid = (lo2 + hi2) // 2
        o, s = block_at(url, mid, size)
        if o is None or s is None or s <= sample:
            lo2 = mid
        else:
            hi2 = mid
    end_blk, _ = block_at(url, hi2, size)
    end = (end_blk if end_blk is not None else size) - 1
    # include one more block of slack so the last sample rows are complete
    end = min(size - 1, end + 200000)
    raw = get(url, start, end)
    sys.stderr.write('range %d-%d (%.1f MB)\n' % (start, end, len(raw) / 1e6))
    # decompress concatenated BGZF members; stop at a truncated trailing member
    data, pos, n_kept = raw, 0, 0
    pref = (sample + '#').encode()
    buf = b''
    with gzip.open(out, 'wb') as fo:
        while pos < len(data):
            d = zlib.decompressobj(31)
            try:
                chunk = d.decompress(data[pos:])
            except zlib.error:
                break
            if not d.eof:
                break
            buf += chunk
            pos = len(data) - len(d.unused_data)
            lines = buf.split(b'\n')
            buf = lines[-1]
            for l in lines[:-1]:
                if l.startswith(pref):
                    fo.write(l + b'\n'); n_kept += 1
    sys.stderr.write('kept %d rows for %s\n' % (n_kept, sample))

if __name__ == '__main__':
    main(*sys.argv[1:4])
