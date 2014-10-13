#!/bin/bash

#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced. 
set -x

###
#code to make script work on both osx and linux. Essentially, it creates a file path to the script directory and saves this path as $0. In detail: 'if the operating system type is darwin (a mac), then use the greadlink function when the readlink function is called. Then use the greadlink function to find the script file named. In doing so, find the path of the script files directory and save as 'scriptdir'. This change is made for Macs because readlink doesn't run properly, but greadlink does. If the OS is not mac (eg. Linux), find the script file using the readlink function and save the path to the script file directory as 'scriptdir.' By using readlink to find where the scripts are, it means if this pipeline is copied onto another computer, the files can still be found.
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
#User defined reference sequence directory.
#This line specifies that the reference directory is located in the script directory in a folder called 'subread_refdir.' The user should create a link for the subread_refdir in the script dir to map to the location of the directory containing their subread indexfiles. These index files can be created by using the 'subread-buildindex' program (refer to RNAseq pipeline user guide) and MUST have the prefix "TAIR10_gen_chrc" (chrc means we included all 7 chromosomes).
refdir=$scriptdir/subread_refdir
#

#Defines the sample that we are working with to the command line as the first token.
sample=$1
sampleFolder=$2

#Specifies the directory that the sample will be opened from. In this case, it is opening a sample folder located in the 'reads_noadapt_trimmed' folder.
sample_dir=${sampleFolder}/$sample

#Defines the output directory to be a folder with the sample name located within the 'align' directory. This will be used in the next step to create an output directory. 
outdir="${sampleFolder}/${sample}"

#Specifies that the sam output file will be placed in the output directory, and have the file name 'sample.sam'
#outsam="${outdir}/${sample}.sam"
outsam="$(ls $sample_dir/*.sam)"

#Specifies that the bam output file will be stored in the output directory, with the file name 'sample.''.bam' has not been added as samtools sort -f currently (29/4/14) has a bug. 
outbam="${outdir}/${sample}" # no .bam, as samtools sort -f has a bug.

#Specifies that the temporary bam output file will be stored in the output directory with the file name 'random.bam.' A temporary bam file has been created due to samtools having a bug with the bam files (Kevins hackery).
tmpbam="${outdir}/${RANDOM}.bam"

echo "samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}"

#Using samtools view to convert the sam file to bam file.
samtools view -S -u $outsam > ${tmpbam}

#Sort the temporary bam file by chromosomal position, and save the sorted file.
samtools sort -m 2G ${tmpbam} $outbam

#Make an index of the sorted bam file
samtools index ${outbam}.bam

#Delete the temporary bam.
rm -v ${outsam} ${tmpbam}
