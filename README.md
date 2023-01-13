## Introduction

vTR (program for computing Variant units in mosaic Tandem Repeats) feeds an input string and an estimated mosaic tandem repeat that partly matches the input, and output a series of variant units.

For example, the upper string in following alignment shows an estimated pattern (AGGGG)3(AAAAGAAAGAGAGGG)2, the lower string is an input string, and the alignment displays 39 matches, 4 mismatches, and 5 indels.

> AGGGG -A-GGGG AGGGG AAAAGAAAGA-GAGGG AAAAGAAAGAGAGGG
>
> AGGGG AAGGGGG A-GGT AAAAGTAAGAGGA-GG AAAAGAAAGAAAGGT

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

