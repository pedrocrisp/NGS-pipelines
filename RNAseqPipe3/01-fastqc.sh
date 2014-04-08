#!/bin/bash

#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced. 
set -x

#If the number of arguments is not one, then echo "USAGE: fastqc.sh SAMPLENAME" and exit the program. Essentially, this code is checking that a sample has been provided and if not, it will exit the program.
if [ $# -ne 1 ]
then
echo "USAGE: fastqc.sh SAMPLENAME"
exit
fi

#Defines the sample that we are working with to the command line as the first token.
sample=$1

#Specifies the directory that the sample will be opened from.
sample_dir=reads/$sample

#Lists all of the files in the particular sample directory with the ending '.fastq.gz' and saves these files as 'fastqs.'
fastqs="$(ls $sample_dir/*.fastq.gz)"

#Creates a new directory called 'reads_qc' and puts the sample names into there. This creates the directory to put the output from fastqc into (next step).
mkdir reads_fastqc/$sample

#Runs the fastq program. This states that the input directory is the fastq files listed in the sample directory (fastqs)and that the output from fastqc goes into the 'reads_fastqc' directory, in a new folder with the sample name. 
fastqc -o reads_fastqc/$sample $fastqs
