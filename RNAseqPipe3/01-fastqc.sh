#!/bin/bash
set -e
set -x
 
if [ $# -ne 1 ]
then
echo "USAGE: fastqc.sh SAMPLENAME"
fi
 
sample=$1
sample_dir=reads/$sample
 
fastqs="$(ls $sample_dir/*.fastq.gz)"
 
mkdir reads_fastqc/$sample
 
fastqc -o reads_fastqc/$sample $fastqs