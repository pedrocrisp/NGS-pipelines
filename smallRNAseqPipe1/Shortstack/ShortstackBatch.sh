#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced. 
set -x

#Defines the sample that we are working with to the command line as the first token.
sample=$1
inputDir=$2
outdir=$3
coverage=$4
reference=$5

#Specifies the directory that the sample will be opened from.
sample_dir="${inputDir}/$sample"

#Find fastQ file
fastQ="$(ls $sample_dir/*q.gz)"
#fastQ="$(ls $sample_dir/*.fastq)"
#Cant get multiple expression to work...
#fastQ="$(ls {${sample_dir}/*q.gz,${sample_dir}/*.fq, ${sample_dir}/*.fastq})"

#Make output directory
outputDir="${outdir}/$sample"
#mkdir $outputDir

#Run Shortstack, default settings, give bowtie some more memory to play with
ShortStack \
--readfile $fastQ \
--sort_mem 4G \
--mincov $coverage \
--genomefile $reference \
--outdir $outputDir

# --readfile a fastq note that must end .fq.gz (not .trimmed.fq.gz - this will break shortstack)
# --sort_mem give samtools more memory (default is 768M.)
# --mincov Deafult: 20 (raw reads) can modify to specify rpm - I have used 5rpm has previously
# --genomefile path to reference genome in .fasta or .fa format. If the bowtie reference is already present this will save time.

#index output bam
samtools index ${outputDir}/${sample}*.bam
