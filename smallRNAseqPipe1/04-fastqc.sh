#!/bin/bash
set -e
set -x
 
if [ $# -ne 1 ]
then
echo "USAGE: fastqc.sh SAMPLENAME"
fi
 
sample=$1
sample_dir=reads_noadapt/$sample
 
fastqs="$(ls $sample_dir/*.fq*)"
 
mkdir reads_noadapt_fastqc/$sample
 
fastqc -o reads_noadapt_fastqc/$sample $fastqs