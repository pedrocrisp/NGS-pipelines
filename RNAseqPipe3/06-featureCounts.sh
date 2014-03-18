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
#reference sequence directory variable - user should create a link called subread_refdir in script dir to the location of the directory containing the subread saf file which must be called "TAIR10_GFF3_genes.saf"
refdir=$scriptdir/subread_refdir
#

sample=$1
sample_dir=align/$sample
outdir="featureCounts/${sample}"
mkdir ${outdir}

featureCounts\
    -F SAF\
    -p\
    -C\
    -s 2\
    -a ${refdir}/TAIR10_GFF3_genes.saf\
    -o "$outdir/${sample}.counts"\
    ${sample_dir}/*.bam
