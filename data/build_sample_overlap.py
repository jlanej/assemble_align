#!/usr/bin/env python3
"""Rebuild the 1kGP cross-resource sample-overlap matrix from public metadata only (no data downloads).

Sources (all HTTPS): NYGC 1kGP 3,202 sequence indices + pedigree, 1kGP phase-3 g1k.ped (3,501), 1kG-ONT S3 bucket listing
(s3://1000g-ont: 500-sample metadata, UltraLong, pod5_data, Illumina_TruPath_WGS, 1KGP_PacBio_WGS, Kinnex, assemblies),
1KG_ONT_VIENNA FTP listing, HPRC release-2 GitHub data tables, HGSVC3 assembly-data table + Hi-C folders,
HGSVC3 2026 Strand-seq pools manifest, HGSVC2 sample list, Platinum Pedigree S3 listing.
Usage: build_sample_overlap.py <out_dir>   (writes sample_overlap.rebuilt.tsv; nothing else is written)
Requires: python>=3.9, pandas, requests. Runtime ~1-6 min (50-70 s measured on 2026-11-11).
Output: the 26 core membership columns of sample_overlap.tsv (2,122 samples present in >=2 resources; 602 complete
NYGC trios). Checked 2026-11-11 against the delivered table: 24 of the 26 columns match cell-for-cell; the assemblies column matches after
stripping parenthetical annotations and '_trio' suffixes, and the rebuilt TruPath boolean (absent from the delivered table)
was checked only against TruPath_cov being non-null. The delivered table adds
23 enrichment columns (per-source coverage, Hi-C, GIAB, T2T contig counts, characterization score, rank, rationale)
that were computed interactively and are not regenerated here.
Notes: three SAMPLE_NAME values in the 2504 sequence index carry trailing spaces (HG02635, HG03025, HG03366) and are
stripped; Kinnex membership = raw data outside 00_quarantine/ or any processed data (HG01271 has only processed data).
"""
import sys, os, re, io, xml.etree.ElementTree as ET
from urllib.parse import urlencode
import requests, pandas as pd
S = requests.Session(); NS = {'s3': 'http://s3.amazonaws.com/doc/2006-03-01/'}
KG = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/'
ONT = 'https://1000g-ont.s3.amazonaws.com'; V = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/'
H3 = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/'; H2 = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/'
HPRC = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/'
PPB = 'https://platinum-pedigree-data.s3.amazonaws.com'
SID = re.compile(r'(HG\d{5}|NA\d{5}|GM\d{5})')
def norm(s): return 'NA' + s[2:] if s.startswith('GM') else s
def get(u):
    r = S.get(u, timeout=180); r.raise_for_status(); return r
def s3_list(bucket, prefix='', delimiter=None):
    tok = None; objs = []; prefs = []
    while True:
        p = {'list-type': '2', 'prefix': prefix}
        if delimiter: p['delimiter'] = delimiter
        if tok: p['continuation-token'] = tok
        root = ET.fromstring(get(bucket + '/?' + urlencode(p)).content)
        objs += [(c.find('s3:Key', NS).text, int(c.find('s3:Size', NS).text)) for c in root.findall('s3:Contents', NS)]
        prefs += [c.find('s3:Prefix', NS).text for c in root.findall('s3:CommonPrefixes', NS)]
        if root.find('s3:IsTruncated', NS).text != 'true': return objs, prefs
        tok = root.find('s3:NextContinuationToken', NS).text
def ls_http(u):
    return [(h, s) for h, d, s in re.findall(r'<a href="([^"?/][^"]*)">[^<]*</a></td><td align="right">([^<]*)</td><td align="right">([^<]*)</td>', get(u).text)]
def walk(u, depth=0, maxdepth=1):
    out = []
    for h, s in ls_http(u):
        if h.endswith('/') and depth < maxdepth: out += walk(u + h, depth + 1, maxdepth)
        else: out.append(u + h)
    return out
def ids(keys): return {norm(m.group(1)) for k in keys for m in [SID.search(k)] if m}

def main(out):
    os.makedirs(out, exist_ok=True)
    def idx(f):
        lines = [l for l in get(KG + f).text.splitlines() if l.strip()]
        hdr = [l for l in lines if l.startswith('#')][-1].lstrip('#').split('\t')
        df = pd.DataFrame([l.split('\t') for l in lines if not l.startswith('#')])
        if df.shape[1] == len(hdr) + 1: df = df.drop(columns=[11])  # 698 index carries one extra empty column
        df.columns = hdr; return df
    ny = pd.concat([idx('1000G_2504_high_coverage.sequence.index').assign(set='2504'), idx('1000G_698_related_high_coverage.sequence.index').assign(set='698')])
    ny['SAMPLE_NAME'] = ny.SAMPLE_NAME.str.strip()  # 3 names in the 2504 index (HG02635, HG03025, HG03366) carry trailing spaces
    ny['cov'] = (pd.to_numeric(ny.BASE_COUNT, errors='coerce') / 3.1e9).round(1)
    nyset = dict(zip(ny.SAMPLE_NAME, ny.set)); nycov = dict(zip(ny.SAMPLE_NAME, ny['cov']))
    pop = pd.read_csv(io.StringIO(get(KG + '20130606_g1k_3202_samples_ped_population.txt').text), sep=r'\s+').set_index('SampleID')
    g1k = pd.read_csv(io.StringIO(get('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped').text), sep='\t').set_index('Individual ID')
    univ = set(pop.index) | set(g1k.index)
    ontm = pd.read_csv(io.StringIO(get(ONT + '/PROCESSED_DATA/1kGP_LRSC_500_ONT_Metadata.tsv').text), sep='\t')
    ont500 = {norm(s): p for s, p in zip(ontm.NHGRI_ID, ontm.ONT_pore)}
    lst = {pre: [k for k, _ in s3_list(ONT, pre)[0]] for pre in ['PROCESSED_DATA/ASSEMBLIES/', 'UltraLong/unaligned_bams/', '1KGP_PacBio_WGS/Raw_Data/', '1KGP_PacBio_WGS/Assembly/', 'Kinnex/kinnex_1KGP/data/raw/', 'Kinnex/kinnex_1KGP/data/processed/']}
    _, tpp = s3_list(ONT, 'Illumina_TruPath_WGS/', '/')
    _, podp = s3_list(ONT, 'pod5_data/', '/')
    pod = {}
    for p in podp:
        m = re.match(r'pod5_data/((?:HG|NA|GM)\d{5})_(R9|R10)(_UL)?/$', p)
        if m: pod.setdefault(norm(m.group(1)), []).append(m.group(2) + ('_UL' if m.group(3) else ''))
    ak = lst['PROCESSED_DATA/ASSEMBLIES/']
    sets = dict(UL=ids(lst['UltraLong/unaligned_bams/']), flye=ids([k for k in ak if '/FLYE/' in k]), r10h=ids([k for k in ak if '/R10_HIFIASM/' in k]),
                shasta=ids([k for k in ak if '/NAPU_PIPELINE_FASTAS/' in k]), pbhifi=ids(lst['1KGP_PacBio_WGS/Raw_Data/']),
                pbtrio=ids([k for k in lst['1KGP_PacBio_WGS/Assembly/'] if '/Trio_Assembly/' in k]), pbasm=ids(lst['1KGP_PacBio_WGS/Assembly/']),
                trupath=ids(tpp), kinnex=ids([k for k in lst['Kinnex/kinnex_1KGP/data/raw/'] if 'quarantine' not in k])
                | ids(lst['Kinnex/kinnex_1KGP/data/processed/']))  # raw outside 00_quarantine, or processed (HG01271: raw quarantined, processed public)
    vienna = {norm(SID.search(h).group(1)) for h, _ in ls_http(V + 'hg38/') if h.endswith('.cram')}
    summ = pd.read_csv(HPRC + 'sequencing_data/data_sample_level_summary_release2_v1.0.csv')
    summ = summ[~summ.sample_id.isin(['CHM13', 'GRCh38'])].set_index('sample_id')
    asm = pd.read_csv(HPRC + 'assemblies_release2_v1.0.index.csv').groupby('sample_id').agg(ph=('phasing', 'first'), me=('assembly_method', 'first'))
    raw = get(H3 + 'release/Assembly_Info/v2.0/supp-table_accessions-assembly-data.2024v2.0.tsv').text.splitlines()
    h = [i for i, l in enumerate(raw) if l.startswith('sample\t')][0]
    a2 = pd.read_csv(io.StringIO('\n'.join(raw[h:])), sep='\t'); a2 = a2[a2.filename == 'all']
    h3 = a2.pivot_table(index='sample', columns='datatype', values='coverage (fold)', aggfunc='first')
    h3['ul'] = a2[a2.datatype == 'ont'].set_index('sample')['coverage >100kbp (fold)']
    hic3 = ids([u.split('/HGSVC3/')[1] for sub in ['working/2021_Hi-C_JAX/', 'working/202505012_USC_HighCovHiC/'] for u in walk(H3 + sub)])
    man = pd.read_csv(io.StringIO(get(H3 + 'working/20260214_Strandseq_Pools/MANIFEST_20260214_Strandseq_Pools.txt').text), sep='\t')
    strand = ids(man.file_path)
    hgsvc2 = ids(re.findall(r'(?:HG|NA|GM)\d{5}', get(H2 + 'HGSVC2_sample_list.txt').text))
    ppids = ids([k for k, _ in s3_list(PPB, '')[0]]) & univ
    pedi = pd.read_csv(io.StringIO(get(KG + 'working/1kGP.3202_samples.pedigree_info.txt').text), sep=r'\s+')  # 3,202-sample pedigree; 602 complete trios once index names are stripped
    trios = {s for s, f, m in zip(pedi.sampleID, pedi.fatherID, pedi.motherID) if f in nyset and m in nyset}
    resources = [set(nyset), set(ont500), vienna, set(summ.index), set(h3.index), sets['pbhifi'], sets['trupath'], strand, sets['kinnex'], ppids, {'NA12878'}]
    def num(x): return None if x is None or pd.isna(x) else float(x)
    rows = []
    for s in sorted(set().union(*resources[1:]) & univ):
        n_res = sum(s in r for r in resources)
        if n_res < 2: continue
        inH, inG = s in summ.index, s in h3.index
        hifi = max([v for v in [num(summ.loc[s, 'coverage_hifi']) if inH else None, num(h3.loc[s, 'hifi']) if inG else None] if v is not None], default=None)
        ul = max([v for v in [num(summ.loc[s, 'coverage_100kb_ont']) if inH else None, num(h3.loc[s, 'ul']) if inG else None] if v is not None], default=None)
        asms = [n for n, c in [('HPRC_R2', inH), ('HGSVC3', inG), ('1kGONT_ShastaHapDup', s in sets['shasta']), ('1kGONT_R10_hifiasm', s in sets['r10h']),
                               ('1kGONT_Flye', s in sets['flye']), ('PacBio1KGP_hifiasm', s in sets['pbasm']), ('PlatinumPed', s in ppids)] if c]
        popn = pop.loc[s, 'Population'] if s in pop.index else g1k.loc[s, 'Population']
        rows.append(dict(sample=s, population=popn, NYGC30x=nyset.get(s, ''), NYGC_raw_cov=nycov.get(s), parents_in_NYGC30x=s in trios,
                         ONT_Miller500=ont500.get(s, ''), ONT_Miller_UL=s in sets['UL'], ONT_pod5=';'.join(sorted(pod.get(s, []))), ONT_Vienna1019_R9=s in vienna,
                         HPRC_R2=inH, HPRC_ONT_chem=summ.loc[s, 'sequencing_chemistry_ont'] if inH else '',
                         HPRC_asm=(asm.loc[s, 'me'] + '/' + asm.loc[s, 'ph']) if s in asm.index else '', HGSVC3=inG, HGSVC3_HiC=s in hic3, StrandSeq=s in strand,
                         PacBio_1KGP_HiFi=s in sets['pbhifi'], PacBio_1KGP_trio_asm=s in sets['pbtrio'], TruPath=s in sets['trupath'], Kinnex_1KGP=s in sets['kinnex'],
                         PlatinumPedigree=s in ppids, HGSVC2_list=s in hgsvc2, HiFi_max_cov=hifi, ONT_100kb_max_cov=ul, assemblies=';'.join(asms),
                         n_assembly_sets=len(asms), n_resources=n_res))
    df = pd.DataFrame(rows).sort_values(['n_resources', 'sample'], ascending=[False, True])
    df.to_csv(os.path.join(out, 'sample_overlap.rebuilt.tsv'), sep='\t', index=False)
    print(f'samples in >=2 resources: {len(df)}; complete NYGC trios: {len(trios)}')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
