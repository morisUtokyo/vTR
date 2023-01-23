#!/bin/bash

# executable modules
uTR=../../uTR/uTR
vTR=../vTR
dir=.
uTR_stat="uTR_stat.txt"

genome="minisatellites_WDR7_ALS"
unit=GTTATGTCCCTTACAGGATTTCACATATCCCTATTCTCAGGCAGGAATAGGGATATGTGAGCTATGATA
TR_file=$dir/$genome".fasta"
TR_uTR=$dir/$genome"_uTR.fasta"
TR_vTR=$dir/$genome"_vTR.fasta"
$uTR -f $TR_file -u $unit -stda -o $TR_uTR 1> $uTR_stat
$vTR -f $TR_uTR > $TR_vTR

genome="minisatellites_ABCA7_Alzheimer"
unit=ACTCCCTCCCCGTGAGCCCCCCCACC
TR_file=$dir/$genome".fasta"
TR_uTR=$dir/$genome"_uTR.fasta"
TR_vTR=$dir/$genome"_vTR.fasta"
$uTR -f $TR_file -u $unit -stda -o $TR_uTR 1> $uTR_stat
$vTR -f $TR_uTR > $TR_vTR

genome="minisatellites_CACNA1C_Bipolar"
unit=AATCACACGACCCTGACCTGACTAGTTTAC
TR_file=$dir/$genome".fasta"
TR_uTR=$dir/$genome"_uTR.fasta"
TR_vTR=$dir/$genome"_vTR.fasta"
$uTR -f $TR_file -u $unit -stda -o $TR_uTR 1> $uTR_stat
$vTR -f $TR_uTR > $TR_vTR

rm $uTR_stat
exit 0
