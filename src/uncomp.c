#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "faidx.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	int i, line_len = 50;
	faidx_t *fai;
	gzFile fp;
	kseq_t *seq;
	if (argc < 3) {
		fprintf(stderr, "Usage: uncomp <ref.fa> <sample.ccomp.fa.gz>\n");
		return 1;
	}
	fai = fai_load(argv[1]);
	fp = gzopen(argv[2], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		char *ref;
		int len, min;
		ref = fai_fetch(fai, seq->name.s, &len);
		min = len < seq->seq.l? len : seq->seq.l;
		printf(">%s", seq->name.s);
		for (i = 0; i < seq->seq.l; ++i) {
			if (i%line_len == 0) putchar('\n');
			if (seq->seq.s[i] == 'Q') seq->seq.s[i] = ref[i];
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
		free(ref);
	}
	kseq_destroy(seq);
	gzclose(fp);
	fai_destroy(fai);
	return 0;
}
