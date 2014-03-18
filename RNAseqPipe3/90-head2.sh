#!/bin/bash
set -e
set -x

sample=$1
sample_dir=reads/$sample

pushd "$sample_dir"
    gunzip *.fastq.gz
    fastqs="$(ls *.fastq)"

    for fq in $fastqs
    do
        fqname="$(basename $fq)"
        outputFile="${fqname%%.*}.small.fastq"
        head -n 40000 $fq > $outputFile
    done

    fastqs="$(ls *.small.fastq)"

    for fq in $fastqs
    do
        fqname="$(basename $fq)"
        outputFile="${fqname%%.*}.small.fastq.gz"
        gzip -c $fq > $outputFile
    done

    rm *.fastq
popd