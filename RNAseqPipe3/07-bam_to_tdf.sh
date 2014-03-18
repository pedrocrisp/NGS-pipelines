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
#reference sequence directory variable - user should create a link called TAIR10_gen_chrc.chrom.sizes in script dir that points to the TAIR10_gen_chrc.genome.sizes file, or put a copy in there (chrc means we included all 7 chromosomes).

chrc_sizes=${scriptdir}/TAIR10_gen_chrc.chrom.sizes


sample=$1
sample_dir=align/$sample
outdir="tdf_for_igv/${sample}"
mkdir ${outdir}

echo "bam to bedgraph"
bedtools genomecov -bg -ibam $sample_dir/*.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes
