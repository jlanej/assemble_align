// Tiny harness: run minimap2 v2.31 mm_sketch() on every sequence of a FASTA file
// and dump (x, y) pairs as raw little-endian uint64 to stdout-specified file.
// x = hash<<8 | span ; y = rid<<32 | lastPos<<1 | strand  (see sketch.c)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "kvec.h"
#include "mmpriv.h"
int main(int argc, char *argv[]) {
	if (argc < 5) { fprintf(stderr, "usage: %s in.fa k w out.bin [hpc]\n", argv[0]); return 1; }
	int k = atoi(argv[2]), w = atoi(argv[3]), hpc = argc > 5 ? atoi(argv[5]) : 0;
	FILE *fi = fopen(argv[1], "r"), *fo = fopen(argv[4], "wb");
	size_t cap = 1<<20, len = 0; char *seq = malloc(cap); char *line = NULL; size_t lc = 0; ssize_t n;
	uint32_t rid = 0; mm128_v p = {0,0,0}; size_t tot = 0;
	#define FLUSH() do { if (len > 0) { p.n = 0; mm_sketch(0, seq, (int)len, w, k, rid, hpc, &p); fwrite(p.a, 16, p.n, fo); tot += p.n; ++rid; len = 0; } } while (0)
	while ((n = getline(&line, &lc, fi)) > 0) {
		if (line[0] == '>') { FLUSH(); continue; }
		while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) --n;
		if (len + n + 1 > cap) { while (len + n + 1 > cap) cap <<= 1; seq = realloc(seq, cap); }
		memcpy(seq + len, line, n); len += n;
	}
	FLUSH();
	fprintf(stderr, "[sketch_harness] k=%d w=%d hpc=%d minimizers=%zu seqs=%u\n", k, w, hpc, tot, rid);
	fclose(fo); fclose(fi); free(seq); free(p.a); return 0;
}
