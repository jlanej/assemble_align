#!/usr/bin/env python3
"""Range-read the first bytes of a remote BAM (e.g. a 1kG TruPath DRAGEN BAM), decode BGZF blocks,
and report @PG/@RG header lines plus auxiliary-tag frequencies over the first N records.
No full download. Stdlib only.
Usage: probe_bam_tags.py --bam-url URL --mb 4 --records 3000 --out tags.json [--scratch DIR]
"""
import argparse, collections, json, os, struct, sys, tempfile, time, urllib.request, zlib

def bgzf_decode(raw):
    buf, off = bytearray(), 0
    while off + 18 <= len(raw) and raw[off:off + 2] == b"\x1f\x8b":
        xlen = struct.unpack("<H", raw[off + 10:off + 12])[0]
        bsize, p = None, off + 12
        while p < off + 12 + xlen:
            si1, si2, slen = raw[p], raw[p + 1], struct.unpack("<H", raw[p + 2:p + 4])[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack("<H", raw[p + 4:p + 6])[0]
            p += 4 + slen
        if bsize is None or off + bsize + 1 > len(raw):
            break
        buf += zlib.decompress(raw[off + 12 + xlen: off + bsize + 1 - 8], -15)
        off += bsize + 1
    return bytes(buf)

SZ = {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}
def aux_tags(b):
    i, out = 0, []
    while i + 3 <= len(b):
        tag, t = b[i:i + 2].decode(errors="replace"), chr(b[i + 2]); i += 3
        if t in "cCA": i += 1
        elif t in "sS": i += 2
        elif t in "iIf": i += 4
        elif t in "ZH": i = b.index(0, i) + 1
        elif t == "B":
            st, cnt = chr(b[i]), struct.unpack("<i", b[i + 1:i + 5])[0]; i += 5 + SZ[st] * cnt; t = "B:" + st
        else: break
        out.append(tag + ":" + t)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-url", required=True); ap.add_argument("--mb", type=int, default=4)
    ap.add_argument("--records", type=int, default=3000); ap.add_argument("--out", required=True)
    ap.add_argument("--scratch", default=None)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory(prefix="scratch.", dir=a.scratch or os.getcwd()) as tmp:
        req = urllib.request.Request(a.bam_url, headers={"Range": f"bytes=0-{a.mb * 1048576 - 1}"})
        raw = urllib.request.urlopen(req, timeout=300).read()
        open(os.path.join(tmp, "head.bin"), "wb").write(raw)
    data = bgzf_decode(raw)
    if data[:4] != b"BAM\x01":
        raise SystemExit("not a BAM")
    lt = struct.unpack("<i", data[4:8])[0]; text = data[8:8 + lt].decode(errors="replace")
    p = 8 + lt; nref = struct.unpack("<i", data[p:p + 4])[0]; p += 4
    for _ in range(nref):
        ln = struct.unpack("<i", data[p:p + 4])[0]; p += 8 + ln
    cnt, n = collections.Counter(), 0
    while p + 4 < len(data) and n < a.records:
        bs = struct.unpack("<i", data[p:p + 4])[0]
        if p + 4 + bs > len(data): break
        rec = data[p + 4:p + 4 + bs]
        l_rn, n_cig, l_seq = rec[8], struct.unpack("<H", rec[12:14])[0], struct.unpack("<i", rec[16:20])[0]
        cnt.update(set(aux_tags(rec[32 + l_rn + 4 * n_cig + (l_seq + 1) // 2 + l_seq:])))
        n += 1; p += 4 + bs
    hl = text.splitlines()
    out = {"bam_url": a.bam_url, "bytes_read": len(raw), "header_lines": len(hl),
           "n_sq": sum(l.startswith("@SQ") for l in hl),
           "pg": [l[:600] for l in hl if l.startswith("@PG")], "rg_example": next((l for l in hl if l.startswith("@RG")), None),
           "records_parsed": n, "records_with_tag": dict(cnt.most_common()),
           "python": sys.version.split()[0], "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(out, open(a.out, "w"), indent=1)
    print(n, "records;", ", ".join(f"{k}={v}" for k, v in cnt.most_common()))

if __name__ == "__main__":
    main()
