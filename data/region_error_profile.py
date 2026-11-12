#!/usr/bin/env python3
"""Per-technology read error profile in one genomic window, using only htslib region slicing over HTTPS.

Reads overlapping REGION are streamed with `samtools view` (no full download), compared base-by-base with
the reference window, and errors are tallied outside +/-MASK bp of the sample's known variants
(taken from a population VCF, e.g. the 1kGP 3,202 phased panel). Homopolymer context = run length of the
reference base at the error position. Output: one TSV row per technology.

Usage:
  region_error_profile.py --region chr20:40000000-40300000 --sample HG01258 \
     --ref-fasta <URL or path to uncompressed FASTA with .fai> \
     --vcf <URL of tabix-indexed VCF containing the sample> \
     --inputs inputs.tsv  (columns: label<TAB>url_of_bam_or_cram) --out out.tsv [--minmapq 20] [--mask 5]
Requires samtools and bcftools (htslib with libcurl) on PATH. Scratch files live in a self-deleting temp dir.
Caveats: errors include residual true variants absent from the VCF, alignment ambiguity, and cell-line
somatic variants; compare technologies on the SAME sample/window, not absolute Q values across studies.
"""
import argparse, subprocess, tempfile, os, re, math, sys, time
from collections import defaultdict

def sh(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(cmd)}\n{r.stderr[-2000:]}")
    return r.stdout

def hp_runs(seq):
    """run length of homopolymer containing each position"""
    n = len(seq); out = [1]*n; i = 0
    while i < n:
        j = i
        while j+1 < n and seq[j+1] == seq[i]: j += 1
        for k in range(i, j+1): out[k] = j-i+1
        i = j+1
    return out

def hpbin(L):
    return '1' if L == 1 else ('2-3' if L <= 3 else ('4-7' if L <= 7 else '>=8'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--region', required=True); ap.add_argument('--sample', required=True)
    ap.add_argument('--ref-fasta', required=True); ap.add_argument('--vcf', required=True)
    ap.add_argument('--inputs', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--minmapq', type=int, default=20); ap.add_argument('--mask', type=int, default=5)
    a = ap.parse_args()
    chrom, rng = a.region.split(':'); start, end = [int(x.replace(',', '')) for x in rng.split('-')]
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as tmp:
        ref_txt = sh(['samtools', 'faidx', a.ref_fasta, f'{chrom}:{start}-{end}'], cwd=tmp)
        ref = ''.join(ref_txt.split('\n')[1:]).upper()
        assert len(ref) == end-start+1, (len(ref), end-start+1)
        hp = hp_runs(ref)
        # sample variants -> masked positions (1-based)
        q = subprocess.run(f"bcftools view -s {a.sample} -r {a.region} {a.vcf} | bcftools view -i 'GT~\"1\"' | bcftools query -f '%POS\\t%REF\\t%ALT\\n'",
                           shell=True, cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if q.returncode != 0: sys.exit(q.stderr[-2000:])
        masked = set(); nvar = 0
        for line in q.stdout.splitlines():
            pos, r_, alt = line.split('\t'); pos = int(pos); nvar += 1
            span = max([len(r_)] + [len(x) for x in alt.split(',')])
            for p in range(pos - a.mask, pos + span + a.mask): masked.add(p)
        inputs = [l.rstrip('\n').split('\t') for l in open(a.inputs) if l.strip() and not l.startswith('#')]
        rows = []
        for label, url in inputs:
            t0 = time.time()
            cmd = ['samtools', 'view', '-F', '0xF04', '-q', str(a.minmapq)]
            if url.endswith('.cram'): cmd += ['-T', a.ref_fasta]
            cmd += [url, a.region]
            p = subprocess.Popen(cmd, cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            nreads = 0; aligned = 0; mm = 0; ins = 0; dele = 0; ins_b = 0; del_b = 0
            by_hp_bases = defaultdict(int); by_hp_indel = defaultdict(int); by_hp_mm = defaultdict(int)
            rlen = []; per_read = []
            for line in p.stdout:
                f = line.split('\t', 11)
                pos = int(f[3]); cig = f[5]; seq = f[9]
                if cig == '*' or seq == '*': continue
                nreads += 1; rlen.append(len(seq))
                rp = pos; qp = 0; r_al = aligned; r_err = mm + ins + dele
                for n_, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
                    n_ = int(n_)
                    if op in 'M=X':
                        for k in range(n_):
                            g = rp + k
                            if start <= g <= end and g not in masked:
                                rb = ref[g-start]
                                if rb == 'N': continue
                                aligned += 1; b = hpbin(hp[g-start]); by_hp_bases[b] += 1
                                if seq[qp+k] != rb and seq[qp+k] != 'N': mm += 1; by_hp_mm[b] += 1
                        rp += n_; qp += n_
                    elif op == 'I':
                        g = rp  # insertion before ref position rp
                        if start <= g <= end and g not in masked and (g-1) not in masked:
                            L = max(hp[g-start], hp[g-1-start] if g-1 >= start else 1)
                            ins += 1; ins_b += n_; by_hp_indel[hpbin(L)] += 1
                        qp += n_
                    elif op in 'DN':
                        if op == 'D':
                            g = rp
                            if start <= g <= end and not any((g+k) in masked for k in range(n_)):
                                dele += 1; del_b += n_; by_hp_indel[hpbin(hp[g-start])] += 1
                        rp += n_
                    elif op == 'S':
                        qp += n_
                ra = aligned - r_al; re_ = (mm + ins + dele) - r_err
                if ra >= 100: per_read.append(re_ / ra)
            p.wait()
            if p.returncode != 0: sys.exit(f"{label}: samtools failed: {p.stderr.read()[-1500:]}")
            err = mm + ins + dele
            rate = err/aligned if aligned else float('nan')
            hp_err = {b: (by_hp_mm[b] + by_hp_indel[b]) / by_hp_bases[b] if by_hp_bases[b] else float('nan') for b in ['1', '2-3', '4-7', '>=8']}
            rlen.sort(); half = sum(rlen)/2; acc = 0; n50 = 0
            for x in reversed(rlen):
                acc += x
                if acc >= half: n50 = x; break
            per_read.sort(); medr = per_read[len(per_read)//2] if per_read else float('nan')
            q20 = sum(1 for x in per_read if x <= 0.01)/len(per_read) if per_read else float('nan')
            q30 = sum(1 for x in per_read if x <= 0.001)/len(per_read) if per_read else float('nan')
            rows.append(dict(label=label, reads=nreads, read_N50=n50, aligned_bases=aligned, mismatch_per_kb=1e3*mm/aligned, ins_events_per_kb=1e3*ins/aligned,
                             del_events_per_kb=1e3*dele/aligned, error_events_per_kb=1e3*rate, phred_event_Q=(-10*math.log10(rate) if rate > 0 else float('inf')),
                             indel_fraction=(ins+dele)/err if err else float('nan'),
                             frac_indels_in_hp4plus=(by_hp_indel['4-7']+by_hp_indel['>=8'])/(ins+dele) if (ins+dele) else float('nan'),
                             err_per_kb_hp1=1e3*hp_err['1'], err_per_kb_hp2_3=1e3*hp_err['2-3'], err_per_kb_hp4_7=1e3*hp_err['4-7'], err_per_kb_hp8plus=1e3*hp_err['>=8'],
                             bases_in_hp8plus=by_hp_bases['>=8'], median_read_Q=(min(60.0, -10*math.log10(medr)) if medr > 0 else 60.0), frac_reads_le_1pct_err=q20, frac_reads_le_0p1pct_err=q30, seconds=round(time.time()-t0, 1)))
            print(f"{label}: reads={nreads} aligned={aligned} err/kb={1e3*rate:.3f} ({time.time()-t0:.0f}s)", file=sys.stderr)
    cols = list(rows[0].keys())
    with open(a.out, 'w') as fo:
        fo.write(f"# region={a.region} sample={a.sample} minmapq={a.minmapq} mask=+/-{a.mask}bp around {nvar} sample variants; masked_positions={len(masked)}\n")
        fo.write('\t'.join(cols) + '\n')
        for r in rows:
            fo.write('\t'.join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + '\n')

if __name__ == '__main__':
    main()
