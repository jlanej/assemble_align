#!/usr/bin/env python3
"""Summarise HPRC release-2 assembly QC (T2T status, gaps, HMM-Flagger HiFi/ONT, NucFlag) per haplotype.
Usage: hprc_r2_qc_baseline.py <repo_raw_base_url> <scratch_dir> <out_dir> [n_io_workers]
Inputs are the public index CSVs in human-pangenomics/hprc_intermediate_assembly/data_tables; every per-haplotype
file they point to (s3://human-pangenomics/...) is fetched over HTTPS into <scratch_dir> and summarised into <out_dir>."""
import sys, os, io, csv, json, hashlib, statistics, collections, concurrent.futures as cf, urllib.request
BASE, SCR, OUT = sys.argv[1].rstrip('/') + '/', sys.argv[2], sys.argv[3]
NW = int(sys.argv[4]) if len(sys.argv) > 4 else 4
IDX = {'t2t': 'data_tables/annotation/chrom_assignment/t2t_sequences_hprc_r2_v1.0.index.csv',
       'gaps': 'data_tables/annotation/chrom_assignment/gaps_hprc_r2_v1.0.index.csv',
       'fl_hifi': 'data_tables/assembly_qc/flagger/flagger_hifi_v0.1.csv',
       'fl_ont': 'data_tables/assembly_qc/flagger/flagger_ont_v0.1.1.csv',
       'nucflag': 'data_tables/assembly_qc/nucflag/nucflag_pre_release_v0.2.index.csv'}
SUMMARY = 'data_tables/sequencing_data/data_sample_level_summary_release2_v1.0.csv'
def get(url, tries=3):
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except Exception as e:
            err = e
    raise err
def https(u): return u.replace('s3://human-pangenomics/', 'https://s3-us-west-2.amazonaws.com/human-pangenomics/')
os.makedirs(SCR, exist_ok=True); os.makedirs(OUT, exist_ok=True)
prov = {'repo_raw_base': BASE, 'index_md5': {}}
index = {}
for k, p in IDX.items():
    b = get(BASE + p); prov['index_md5'][p] = hashlib.md5(b).hexdigest()
    index[k] = list(csv.DictReader(io.StringIO(b.decode())))
sb = get(BASE + SUMMARY); prov['index_md5'][SUMMARY] = hashlib.md5(sb).hexdigest()
summ = {r['sample_id']: r for r in csv.DictReader(io.StringIO(sb.decode()))}
jobs = []
for k, rows in index.items():
    for r in rows:
        jobs.append((k, r['assembly_name'], https(r['location'])))
def fetch(job):
    k, asm, url = job
    fn = os.path.join(SCR, f'{asm}.{k}')
    if not os.path.exists(fn):
        try:
            data = get(url)
        except Exception as e:
            return (k, asm, None, str(e)[:120])
        open(fn, 'wb').write(data)
    return (k, asm, fn, None)
files = collections.defaultdict(dict); fails = []
with cf.ThreadPoolExecutor(NW) as ex:
    for k, asm, fn, err in ex.map(fetch, jobs):
        if fn: files[asm][k] = fn
        else: fails.append((k, asm, err))
def read_bed(fn, skip_classes=()):
    d = collections.defaultdict(list)
    for line in open(fn):
        if not line.strip() or line.startswith(('track', '#')): continue
        f = line.rstrip('\n').split('\t')
        c = f[3] if len(f) > 3 else '.'
        if c in skip_classes: continue
        d[f[0]].append((int(f[1]), int(f[2]), c))
    return d
def merge(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
def total(m): return sum(e - s for s, e in m)
def inter(a, b):
    i = j = 0; t = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e > s: t += e - s
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return t
rows = []; chr_t2t = collections.Counter(); chr_n = collections.Counter(); nuc_classes = collections.Counter()
for asm, fs in sorted(files.items()):
    sample = asm.split('_')[0]; hap = asm.split('_')[1]
    r = {'assembly_name': asm, 'sample_id': sample, 'hap': hap}
    if 't2t' in fs:
        nc = ns = 0
        for rec in csv.DictReader(open(fs['t2t']), delimiter='\t'):
            if rec['level'] == 'contig': nc += 1
            else: ns += 1
            chr_t2t[(rec['chr_name'], rec['level'])] += 1
        r['t2t_contigs'] = nc; r['t2t_scaffolds_with_gaps'] = ns; r['t2t_total'] = nc + ns
    if 'gaps' in fs:
        g = [l.split('\t') for l in open(fs['gaps']) if l.strip()]
        r['n_gaps'] = len(g); r['gap_bp'] = sum(int(x[2]) - int(x[1]) for x in g)
    merged = {}
    for k in ('fl_hifi', 'fl_ont'):
        if k in fs:
            bed = read_bed(fs[k], skip_classes=('NNN', 'Hap'))
            for cls in ('Err', 'Dup', 'Col'):
                r[f'{k}_{cls}_bp'] = sum(total(merge([(s, e) for s, e, c in v if c == cls])) for v in bed.values())
            merged[k] = {ctg: merge([(s, e) for s, e, c in v]) for ctg, v in bed.items()}
            r[f'{k}_flagged_bp'] = sum(total(m) for m in merged[k].values())
    if 'fl_hifi' in merged and 'fl_ont' in merged:
        both = sum(inter(merged['fl_hifi'].get(c, []), merged['fl_ont'].get(c, [])) for c in set(merged['fl_hifi']) | set(merged['fl_ont']))
        r['fl_hifi_and_ont_bp'] = both
    if 'nucflag' in fs:
        bed = read_bed(fs['nucflag'])
        for v in bed.values():
            for s, e, c in v: nuc_classes[c] += e - s
        nonhet = {ctg: merge([(s, e) for s, e, c in v if c not in ('HET', 'GOOD', 'NNN', 'GAP')]) for ctg, v in bed.items()}
        r['nucflag_nonhet_bp'] = sum(total(m) for m in nonhet.values())
        r['nucflag_misjoin_bp'] = sum(total(merge([(s, e) for s, e, c in v if c == 'MISJOIN'])) for v in bed.values())
        if 'fl_hifi' in merged:
            un = {c: merge(merged['fl_hifi'].get(c, []) + merged.get('fl_ont', {}).get(c, [])) for c in set(merged['fl_hifi']) | set(merged.get('fl_ont', {}))}
            r['nucflag_nonhet_and_flagger_any_bp'] = sum(inter(nonhet.get(c, []), un.get(c, [])) for c in set(nonhet) | set(un))
    s = summ.get(sample, {})
    for col in ('coverage_hifi', 'coverage_ont', 'coverage_100kb_ont', 'sequencing_chemistry_ont', 'coverage_hic', 'read_N50_ont', 'n50_hifi'):
        r[col] = s.get(col, '')
    rows.append(r)
cols = sorted({k for r in rows for k in r}, key=lambda x: (x not in ('assembly_name', 'sample_id', 'hap'), x))
with open(os.path.join(OUT, 'hprc_r2_qc_per_haplotype.tsv'), 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader(); [w.writerow(r) for r in rows]
with open(os.path.join(OUT, 'hprc_r2_t2t_by_chrom.tsv'), 'w') as fh:
    fh.write('chr\tlevel\tn_haplotypes\n')
    for (c, l), n in sorted(chr_t2t.items()): fh.write(f'{c}\t{l}\t{n}\n')
prov.update({'n_assemblies_summarised': len(rows), 'n_fetch_failures': len(fails), 'fetch_failures': fails[:50],
             'nucflag_class_bp_total': dict(nuc_classes)})
json.dump(prov, open(os.path.join(OUT, 'provenance.json'), 'w'), indent=1)
print('assemblies', len(rows), 'failures', len(fails))
