## Real data

The following minisatellite DNA sequences (of units >20bp in size) that are expanded in case samples of three brain disorders, ALS, Alzheimer's disease, and Bipolar disorder. 

- WDR7(control,ALS) hg38_dna range=chr18:57024495-57024955 unit=GTTATGTCCCTTACAGGATTTCACATATCCCTATTCTCAGGCAGGAATAGGGATATGTGAGCTATGATA

- ABCA7(control,Alzheimer) hg38_dna range=chr19:1049431-1050033 unit=ACTCCCTCCCCGTGAGCCCCCCCACC

- CACNA1C(control,Bipolar) hg38_dna range=chr12:2255791-2256090 unit=AATCACACGACCCTGACCTGACTAGTTTAC

As DNA sequences from case samples are not availabe for privacy protection reasons, we put minisatellite sequences collected from the reference genome hg38 into the following files:

    minisatellites_WDR7_ALS.fasta
    minisatellites_ABCA7_Alzheimer.fasta
    minisatellites_CACNA1C_Bipolar.fasta

In the above files, annotations describe genes (e.g., WDR7) associated with individual disorders, the ranges of position (e.g., chr18:57024495-57024955) in the reference genome hg38, and the representative units that are moderately expanded even in the reference genome.
The representative units are of length >25 nt and are not identified by uTR that searchs for units of length at most 20 nt. The above long units are detected using our software program named mTR, which is available at:

    https://github.com/morisUtokyo/mTR
    
To calculate variant unis in mosaic tandem repeats, we use the script program "test.sh," which assume that uTR is in the directory ../../uTR/uTR and vTR is in ../vTR/vTR.

"test.sh" first applies uTR to decompose the DNA sequences in each of the above fasta files into mosaic tandem repeats, and put the result into: 

    minisatellites_WDR7_ALS_uTR.fasta
    minisatellites_ABCA7_Alzheimer_uTR.fasta
    minisatellites_CACNA1C_Bipolar_uTR.fasta

Each of the above files stores minisatellite patterns:

    <GTTATGTCCCTTACAGGATTTCACATATCCCTATTCTCAGGCAGGAATAGGGATATGTGAGCTATGATA>6 for WDR7(control,ALS)
    <ACTCCCTCCCCGTGAGCCCCCCCACC>23 for ABCA7(control,Alzheimer)
    <AATCACACGACCCTGACCTGACTAGTTTAC>10 for CACNA1C(control,Bipolar)

Afterwards, using vTR, "test.sh" calcualtes variant units and put them into:

    minisatellites_WDR7_ALS_vTR.fasta
    minisatellites_ABCA7_Alzheimer_vTR.fasta
    minisatellites_CACNA1C_Bipolar_vTR.fasta
