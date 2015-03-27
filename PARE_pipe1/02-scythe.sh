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

###User to ensure correct file path#

#Scythe relies upon knowing the reference directory for the adaptors used, so that it can recognise these sequences. 
#This code is defining the file path to the adaptor sequences, located within the script directory. 
adapterfile="$scriptdir/truseq_adapters.fasta"

#Defines the sample that we are working with to the command line as the first token.
sample=$1

#Specifies the directory that the sample will be opened from.
sample_dir=reads/$sample

#List all files ending with '.fastq.gz' that are located within the specified sample directory and save these as the variable 'fastqs.'
fastqs="$(ls $sample_dir/*.fastq.gz)"

#Creates a new directory called 'reads_noadapt' and within this a folder for the sample. This creates the directory to put the output from scythe into (next step).
outdir="reads_noadapt_scythe/$sample"
mkdir $outdir

#This command runs scythe. It says 'for the fastqs listed within the sample directory do the following:
#1) keep the file names but remove the file extensions
#2) Store the output files in the specified sample folder within the reads_noadapt directory. Store the file name as 'samplename.noadapt.fq.gz' 
#3) Run the scythe function
#4) -p: Change the prior from the default (0.3) to 0.1
#5) -a: The adaptor sequence is located in the script directory, and is called 'truseq_adapters.fasta'
#6) run scythe on this given fastq file.
#7) move the output of scythe to the output file (reads_noadapt/given sample name)

for fq in $fastqs
do
fqname="$(basename $fq)"
outputFile="$outdir/${fqname%%.*}.noadapt.fq"
scythe \
-p 0.1 \
-a $adapterfile \
$fq \
>$outputFile
done

#gzip fastq (differs from RNAseq pipe because fastqs are not passed to quality trimmer)

gzip "$outdir/${fqname%%.*}.noadapt.fq"
