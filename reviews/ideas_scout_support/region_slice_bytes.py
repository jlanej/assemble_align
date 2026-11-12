#!/usr/bin/env python3
"""Estimate how many bytes of a remote CRAM must be fetched to retrieve all reads that
overlap a set of regions, using only the CRAM's .crai index (no reads are downloaded).

Contig names for .crai sequence ids are read from the CRAM header (first container), so
the mapping is not assumed. Region sets are BED files (plain or gzipped, URL or path) or
UCSC genomicSuperDups.txt(.gz) tables. Overlapping intervals are merged per set and can be
padded. Bytes = sum of .crai slice sizes for slices that overlap any merged interval
(container headers are ignored, so the figure is a slight underestimate). Each region set is
reported for two scopes: 'primary' (chr1-22, X, Y, M by default) and 'all' contigs in the
CRAM header (adds alt, random, unplaced and HLA contigs present in the region table).

Example:
  region_slice_bytes.py --cram-url https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram \
     --ucsc-sd SD=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
     --bed MRJD=https://1000g-ont.s3.amazonaws.com/Illumina_TruPath_WGS/HG00097/MRJD_regions_38.bed.gz \
     --pad 0 --pad 1000 --out result.json --scratch "$TMPDIR"
Stdlib only (Python >= 3.8).
"""
import argparse, gzip, json, os, re, shutil, sys, tempfile, time, urllib.request, zlib

def fetch(src, dest):
    if re.match(r"https?://", src):
        with urllib.request.urlopen(src, timeout=300) as r, open(dest, "wb") as fh:
            shutil.copyfileobj(r, fh)
    else:
        shutil.copy(src, dest)
    return dest

def read_text(path):
    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode(errors="replace").splitlines()

def cram_contigs(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-1048575"})
    raw = urllib.request.urlopen(req, timeout=120).read()
    if raw[:4] != b"CRAM":
        raise SystemExit("not a CRAM file: " + url)
    i = raw.find(b"@HD")
    txt = raw[i:].decode(errors="replace") if i >= 0 else None
    if txt is None:
        j = raw.find(b"\x1f\x8b")
        txt = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw[j:]).decode(errors="replace")
    return re.findall(r"@SQ\tSN:([^\t\n]+)", txt), "%d.%d" % (raw[4], raw[5])

def merge(iv):
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cram-url", required=True)
    ap.add_argument("--md5", default=None, help="expected CRAM md5 (recorded only)")
    ap.add_argument("--bed", action="append", default=[], help="LABEL=URL_OR_PATH (BED, optionally gzipped)")
    ap.add_argument("--ucsc-sd", action="append", default=[], help="LABEL=URL_OR_PATH (UCSC genomicSuperDups table)")
    ap.add_argument("--pad", action="append", type=int, default=[], help="padding in bp (repeatable); default 0")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scratch", default=None, help="parent directory for a self-cleaning temp dir")
    ap.add_argument("--primary-regex", default=r"^chr([0-9]+|X|Y|M)$",
                    help="contigs counted in the 'primary' scope; the 'all' scope uses every contig in the CRAM header")
    a = ap.parse_args()
    pads = a.pad or [0]
    with tempfile.TemporaryDirectory(prefix="scratch.", dir=a.scratch or os.getcwd()) as tmp:
        head = urllib.request.urlopen(urllib.request.Request(a.cram_url, method="HEAD"), timeout=120)
        cram_bytes = int(head.headers["Content-Length"])
        names, version = cram_contigs(a.cram_url)
        crai = read_text(fetch(a.cram_url + ".crai", os.path.join(tmp, "x.crai")))
        by = {}
        for ln in crai:
            if not ln.strip():
                continue
            sid, st, sp, co, so, ss = (int(x) for x in ln.split("\t")[:6])
            if sid < 0:
                continue
            by.setdefault(names[sid], []).append((st, st + sp, co, so, ss))
        sets = {}
        for spec in a.bed:
            lab, src = spec.split("=", 1)
            rows = [l.split("\t") for l in read_text(fetch(src, os.path.join(tmp, lab + ".bed")))
                    if l and not l.startswith(("#", "track", "browser"))]
            sets[lab] = {"source": src, "iv": [(r[0], int(r[1]), int(r[2])) for r in rows]}
        for spec in a.ucsc_sd:
            lab, src = spec.split("=", 1)
            rows = [l.split("\t") for l in read_text(fetch(src, os.path.join(tmp, lab + ".txt"))) if l]
            sets[lab] = {"source": src, "iv": [(r[1], int(r[2]), int(r[3])) for r in rows]}
        res = {}
        prim = re.compile(a.primary_regex)
        for lab, d in sets.items():
          for scope in ("primary", "all"):
            for pad in pads:
                m = merge([(c, max(0, s - pad), e + pad) for c, s, e in d["iv"]])
                m = [x for x in m if x[0] in by and (scope == "all" or prim.match(x[0]))]
                sel = set()
                for c, s, e in m:
                    for (x0, x1, co, so, ss) in by[c]:
                        if x0 <= e and x1 >= s:
                            sel.add((co, so, ss))
                b = sum(x[2] for x in sel)
                res[f"{lab}_{scope}_pad{pad}"] = {"merged_intervals": len(m), "merged_bp": sum(e - s for _, s, e in m),
                                          "slices": len(sel), "bytes": b, "fraction_of_cram": b / cram_bytes}
    out = {"cram_url": a.cram_url, "cram_md5_expected": a.md5, "cram_version": version, "cram_bytes": cram_bytes,
           "crai_entries": len(crai), "n_contigs_in_header": len(names), "region_sets": {k: v["source"] for k, v in sets.items()},
           "results": res, "python": sys.version.split()[0], "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "method": "sum of .crai slice sizes for slices overlapping merged (padded) intervals; container headers ignored"}
    json.dump(out, open(a.out, "w"), indent=1)
    for k, v in res.items():
        print(f"{k}\t{v['merged_bp']}\t{v['bytes']}\t{100 * v['fraction_of_cram']:.2f}%")

if __name__ == "__main__":
    main()
