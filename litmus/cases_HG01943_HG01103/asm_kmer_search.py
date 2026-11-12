import sys, gzip, urllib.request, ahocorasick, json, io, time
"""Stream FASTA(s) (URL or path, gzip/bgzip or plain) and report exact hits of k-mer patterns (both strands).
usage: asm_kmer_search.py patterns.tsv(name<TAB>seq) out.json fasta1 [fasta2 ...]"""
pat_f, out_f, fas = sys.argv[1], sys.argv[2], sys.argv[3:]
comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
A = ahocorasick.Automaton()
for line in open(pat_f):
    n, s = line.split()[:2]; s = s.upper()
    A.add_word(s, (n, '+', len(s))); A.add_word(s.translate(comp)[::-1], (n, '-', len(s)))
A.make_automaton()
hits = {}
def scan(name, seq, src):
    for end, (n, strand, k) in A.iter(seq):
        hits.setdefault(n, []).append(dict(source=src, contig=name, pos=end - k + 2, strand=strand))
for fa in fas:
    t0 = time.time()
    fh = urllib.request.urlopen(fa, timeout=600) if fa.startswith('http') else open(fa, 'rb')
    if fa.endswith('.gz'): fh = gzip.GzipFile(fileobj=fh)
    txt = io.TextIOWrapper(fh, encoding='ascii', errors='replace')
    name, buf, nseq, nbp = None, [], 0, 0
    src = fa.split('/')[-1]
    for line in txt:
        if line.startswith('>'):
            if name is not None:
                s = ''.join(buf).upper(); scan(name, s, src); nseq += 1; nbp += len(s)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        s = ''.join(buf).upper(); scan(name, s, src); nseq += 1; nbp += len(s)
    print(f'{src}: {nseq} sequences, {nbp} bp scanned in {time.time()-t0:.0f}s', flush=True)
json.dump(hits, open(out_f, 'w'), indent=1)
summary = {}
for line in open(pat_f):
    n = line.split()[0]; summary[n] = len(hits.get(n, []))
print(json.dumps(summary))
