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
coverage=$5

sample_dir=$alignFolder/$sample
outFolder="${outdir}/${sample}"
mkdir ${outFolder}

####### non-stranded module


####### stranded module
if [ "$coverage" == "5prime" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.plus5.bed > $outFolder/$sample.plus.real.bed
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.minus5.bed > $outFolder/$sample.minus.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.plus.real.bed > $outFolder/$sample.plus_sorted.bed
sort -k1,1 -k2,2n $outFolder/$sample.minus.real.bed > $outFolder/$sample.minus_sorted.bed
#remove intermediate file
rm -rv $outFolder/*real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.plus_sorted.bed -b $reference > $outFolder/$sample.plus.dist.bed
closestBed -D "ref" -a $outFolder/$sample.minus_sorted.bed -b $reference > $outFolder/$sample.minus.dist.bed
#remove intermediate file
rm -rv $outFolder/*_sorted.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.plus.dist.bed > $outFolder/$sample.plus.dist.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.minus.dist.bed > $outFolder/$sample.minus.dist.1k.bed

elif [ "$coverage" == "full" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.plus.bed > $outFolder/$sample.plus.real.bed
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.minus.bed > $outFolder/$sample.minus.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.plus.real.bed > $outFolder/$sample.plus_sorted.bed
sort -k1,1 -k2,2n $outFolder/$sample.minus.real.bed > $outFolder/$sample.minus_sorted.bed
#remove intermediate file
rm -rv $outFolder/*real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.plus_sorted.bed -b $reference > $outFolder/$sample.plus.dist.bed
closestBed -D "ref" -a $outFolder/$sample.minus_sorted.bed -b $reference > $outFolder/$sample.minus.dist.bed
#remove intermediate file
rm -rv $outFolder/*_sorted.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.plus.dist.bed > $outFolder/$sample.plus.dist.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.minus.dist.bed > $outFolder/$sample.minus.dist.1k.bed

elif [ "$coverage" == "nonstranded" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $sample_dir/$sample.bed > $outFolder/$sample.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.real.bed > $outFolder/$sample.sorted.bed

#remove intermediate file
rm -rv $outFolder/*real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.sorted.bed -b $reference > $outFolder/$sample.dist.bed

#remove intermediate file
rm -rv $outFolder/*sorted.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.dist.bed > $outFolder/$sample.dist.1k.bed

else
echo "ERROR: coverage type has not been specificed select 5prime or full"
exit 1
fi

# Then remove any files absolutely required from making the plot, Just keep the dist.1k.bed
rm -rv $outFolder/*dist.bed
