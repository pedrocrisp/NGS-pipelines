#!/bin/bash
set -e
set -x

sample=$1
sample_dir=reads_scythe/$sample
 
fastqs="$(ls $sample_dir/*.fq.gz)"
 
mkdir reads_scythe_seqtk/$sample

for fq in $fastqs
do
fqname="$(basename $fq)"
outputFile="reads_scythe_seqtk/$sample/${fqname%%.*}.trimmed.fq"
seqtk trimfq \
-l 1 \
$fq \
>$outputFile
done

