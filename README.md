## Introduction

vTR (program for computing Variant units in mosaic Tandem Repeats) feeds an input string and an estimated mosaic tandem repeat that partly matches the input, and output a series of variant units.

For example, the upper string in following alignment shows an estimated pattern (AGGGG)3(AAAAGAAAGAGAGGG)2, the lower string is an input string, and the alignment displays 39 matches, 4 mismatches, and 5 indels.

> AGGGG -A-GGGG AGGGG AAAAGAAAGA-GAGGG AAAAGAAAGAGAGGG
> 
> AGGGG AAGGGGG A-GGT AAAAGTAAGAGGA-GG AAAAGAAAGAAAGGT

#include <string.h>
#include <stdio.h>
#include "ksw2.h"

void align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	putchar('\n');
	free(ez.cigar); free(ts); free(qs);
}

int main(int argc, char *argv[])
{
	align("ATAGCTAGCTAGCAT", "AGCTAcCGCAT", 1, -2, 2, 1);
	return 0;
}



Although the estimated pattern, (AGGGG)3(AAAAGAAAGAGAGGG)2, is concise and is easier to understand the structure of the mosaic tandem repeat, to see what substitutions and indels are present in the lower string precisely, it would be useful to examine a series of unit variants in the string:

> AGGGG AAGGGGG AGGT AAAAGTAAGAGGAGG AAAAGAAAGAAAGGT

For this purpose, vTR is developed to output the above decomposition from the estimated pattern and the input string. vTR feeds a fasta file of the form

> \> #Pat (AGGGG)3(AAAAGAAAGAGAGGG)2
>
> AGGGGAAGGGGGAGGTAAAAGTAAGAGGAGGAAAAGAAAGAAAGGT

where the pattern following #Pat shows the original pattern, and vYT outputs

> \> #Len 46 #Err 0.188 #Pat <AGGGG>3<AAAAGAAAGAGAGGG>2 #PrecisePat <AGGGG>1<AAGGGGG>1<AGGT>1<AAAAGTAAGAGGAGG>1<AAAAGAAAGAAAGGT>1
>
> AGGGGAAGGGGGAGGTAAAAGTAAGAGGAGGAAAAGAAAGAAAGGT

where "#Len 46" means the input is of length 46, and "#Err 0.188" shows the error rate between the pattern and the input string. 

The error rate is defined as the sum of mismatches, insertions, and deletions (denoted by X) devided by the sum of the number of matches and X (denoted by Y), namely X/Y. In the running example, X = 9(=4+5), Y = 48, and X/Y = 0.188 = 9/48.

A sample code with data is available at test/test_sample.sh.

vTR use KSW2, a library to align a pair of reads that implements a global alignment dynamic programming algorithm. You can obtain a copy of the KSW2 program from:

https://github.com/lh3/ksw2

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100.

For your convenience, a copy of the KSW2 program is placed on this github.

## Usage

vTR -f fasta_file] 
* -f: Input a fasta file, say sample.fasta

