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

###
#agrs
sample=$1
alignFolder=$2
reference=$3
outdir=$4
sample_dir=$alignFolder/$sample
outFolder="${outdir}/${sample}"
mkdir ${outFolder}

#make a real bed file by adding chromosome edn co-ordinates (same as start because nt resolution)
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.plus.bed > $outFolder/$sample.plus.real.bed
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.minus.bed > $outFolder/$sample.minus.real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.plus.real.bed -b $reference > $outFolder/$sample.plus.dist.bed
closestBed -D "ref" -a $outFolder/$sample.minus.real.bed -b $reference > $outFolder/$sample.minus.dist.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.plus.dist.bed > $outFolder/$sample.plus.dist.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.minus.dist.bed > $outFolder/$sample.minus.dist.1k.bed





