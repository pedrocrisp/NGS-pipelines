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

#Defines the sample that we are working with to the command line as the first token.
sample=$1

#Specifies the directory that the sample will be opened from. In this case, it is opening a sample folder located in the 'reads_noadapt_trimmed' folder.
sample_dir=align/$sample

#Defines the output directory to be a folder with the sample name located within the 'align' directory. This will be used in the next step to create an output directory. 
outdir="align_rmdup/${sample}"

#Creates an output directory to put the returned files to go into once subread has been run on the sample. In this case, the output from subread for the sample should go into a folder containing the sample's name, located within the 'align' directory.
mkdir ${outdir}

#List the file ending with '.bam' that are located within the specified sample directory and save as the variable 'bamname.'

bamname="$(ls $sample_dir/*.bam)"

#make output file name
outbamname="$(basename $bamname)"
outputFile="${outdir}/${outbamname%%.*}_rmdup.bam"

samtools rmdup $bamname $outputFile
