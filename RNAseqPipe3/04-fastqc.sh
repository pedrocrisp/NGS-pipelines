#!/bin/bash
set -e
set -x
 
if [ $# -ne 1 ]
then
echo "USAGE: fastqc.sh SAMPLENAME"
fi
 
sample=$1
sample_dir=reads_scythe_seqtk/$sample
 
fastqs="$(ls $sample_dir/*.fq)"
 
mkdir reads_scythe_seqtk_fastqc/$sample
 
fastqc -o reads_scythe_seqtk_fastqc/$sample $fastqs