#!/bin/bash
set -e
set -x

sample=$1
inputDir=reads_noadapt
sample_dir="${inputDir}/$sample"
outDir=reads_noadapt_trimmed

forward_fq="$(ls $sample_dir/*R1.noadapt.fq.gz)"
reverse_fq="$(ls $sample_dir/*R2.noadapt.fq.gz)"

mkdir "${outDir}/$sample"


forward_fqname="$(basename $forward_fq)"
forward_fq_outputFile="${outDir}/$sample/${forward_fqname%%.*}.trimmed.fq"

reverse_fqname="$(basename $reverse_fq)"
reverse_fq_outputFile="${outDir}/$sample/${reverse_fqname%%.*}.trimmed.fq"


sickle pe \
-f $forward_fq \
-r $reverse_fq \
-t sanger \
-o $forward_fq_outputFile \
-p $reverse_fq_outputFile \
-s "${outDir}/$sample/${sample}.trimmed.singles.fq" \
-q 20 \
-l 20


