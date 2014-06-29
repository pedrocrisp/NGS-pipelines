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

#split F and R read
samtools view -F 64  -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam
samtools view -f 64 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.reverse.bam

echo "bam to bedgraph"
#non-stratded bedgraph
bedtools genomecov -bg -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#forward reads stranded
bedtools genomecov -strand + -ibam $sample_dir/*forward.bam -g $chrc_sizes -d > $outdir/${sample}_forward_plus.bed

bedtools genomecov -strand - -ibam $sample_dir/*forward.bam -g $chrc_sizes -d > $outdir/${sample}_forward_minus.bed
#-scale -1

#reverse reads stranded
bedtools genomecov -strand + -ibam $sample_dir/*reverse.bam -g $chrc_sizes -d > $outdir/${sample}_reverse_plus.bed
#-scale -1

bedtools genomecov -strand - -ibam $sample_dir/*reverse.bam -g $chrc_sizes -d > $outdir/${sample}_reverse_minus.bed

#merge F read on + strand and R read on - strand

cat $outdir/*forward_plus.bed $outdir/*reverse_minus.bed | mergeBed -i stdin > $outdir/${sample}.plus.bed

#merge F read on - strand and R read on + strand

cat $outdir/*reverse_plus.bed $outdir/*forward_minus.bed | mergeBed -i stdin > $outdir/${sample}.minus.bed

#make bedgraphs
bedtools genomecov -bg -i $outdir/${sample}.plus.bed -g $chrc_sizes > $outdir/${sample}.plus.bg

bedtools genomecov -bg -scale -1 -i $outdir/${sample}.minus.bed -g $chrc_sizes > $outdir/${sample}.minus.bg

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWig
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

sort -k1,1 -k2,2n $outdir/${sample}.minus.bg > $outdir/${sample}.sorted.minus.bg

sort -k1,1 -k2,2n $outdir/${sample}.plus.bg > $outdir/${sample}.sorted.plus.bg

bedGraphToBigWig $outdir/*.sorted.plus.bg $chrc_sizes $outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/*.sorted.minus.bg $chrc_sizes $outdir/$sample.minus.bigWig
