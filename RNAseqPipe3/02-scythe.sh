#!/bin/bash
set -e
set -x

###
#code to make it work on osx and linux
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

adapterfile="$scriptdir/truseq_adapters.fasta"

sample=$1
sample_dir=reads/$sample
 
fastqs="$(ls $sample_dir/*.fastq.gz)"
 
mkdir reads_scythe/$sample

for fq in $fastqs
do
fqname="$(basename $fq)"
outputFile="reads_scythe/$sample/${fqname%%.*}.noadapt.fq.gz"
scythe \
-p 0.1 \
-a $adapterfile \
$fq \
>$outputFile
done

# -p set prior to 0.1 default 0.3
# -a adapter file
# -o output file
#input