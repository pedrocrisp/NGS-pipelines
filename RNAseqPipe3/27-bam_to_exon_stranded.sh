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


####################
##code to create exon reference beds:
##get exons from GFF
#
#grep exon TAIR10_GFF3_genes.gff >TAIR10_GFF3_exons.gff
#
##make bed
#gff2bed <TAIR10_GFF3_exons.gff >TAIR10_GFF3_exons.bed
#
##fix bed
#cut -f 10 TAIR10_GFF3_exons.bed | cut -f 2 -d '=' |paste TAIR10_GFF3_exons.bed -  > temp_TAIR10_GFF3_exons.bed
#
#cut -f 1-3,5- temp_TAIR10_GFF3_exons.bed > temp2_TAIR10_GFF3_exons.bed
#
#awk -F'\t' -v OFS="\t" '{$4=$NF OFS $4;$NF=""}7' temp2_TAIR10_GFF3_exons.bed > TAIR10_exons.bed
#
#grep + TAIR10_exons.bed > TAIR10_exons_plusStrand.bed
#grep - TAIR10_exons.bed > TAIR10_exons_minusStrand.bed
#
####################
#
exonsPlusBed=${scriptdir}/TAIR10_exons_plusStrand.bed
exonsMinusBed=${scriptdir}/TAIR10_exons_minusStrand.bed

sample=$1
sample_dir=align/$sample
outdir="exon_beds/${sample}"
mkdir ${outdir}

######################
####The bams must be split into the "plus strand reads and minus strands reads in seperate files so that bedtools can be run non-stranded in order to create strandedness because bedtools doent know about PE reads... below code creates the necessary bam, this sould be done in step 07-
##split F and R reads into plus and minus strand taking into account PE
##http://seqanswers.com/forums/showthread.php?t=29399
#
##R1 forward strand
#samtools view -f 99 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1F.bam
#
##R2 reverse strand
#samtools view -f 147 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2R.bam
#
#samtools merge -f $sample_dir/${sample}.forward.bam $sample_dir/${sample}.R1F.bam $sample_dir/${sample}.R2R.bam
#
##R1 reverse strand
#samtools view -f 83 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1R.bam
#
##R2 forward strand
#samtools view -f 163 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2F.bam
#
#samtools merge -f $sample_dir/${sample}.reverse.bam $sample_dir/${sample}.R1R.bam $sample_dir/${sample}.R2F.bam
#
#rm $sample_dir/${sample}*.R*.bam
#
#####################

echo "bam to exon-bedgraph"
##exon converage beds
#plus strand (dont include -s flag!)
coverageBed -abam $sample_dir/*reverse.bam -b $exonsPlusBed -d > $sample.plus.exons.bed

#minus strand
coverageBed -abam $sample_dir/*forward.bam -b $exonsMinusBed -d > $sample.minus.exons.bed


