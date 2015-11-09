#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced. 
set -x

#Defines the sample that we are working with to the command line as the first token.
sample=$1
inputDir=$2
outdir=$3

#Specifies the directory that the sample will be opened from.
sample_dir="${inputDir}/$sample"

#Find fastQ file
fastQ="$(ls $sample_dir/*.fq.gz)"

#Make output directory
outputDir="${outDir}/$sample"
mkdir $outputDir

#Run Shortstack, default settings, give bowtie some more memory to play with
ShortStack \
--readfile $fastQ \
--sort_mem 4G \
--genomefile ~/ws/refseqs/TAIR10/chromosomes/TAIR10.fa \
--outdir $outputDir




