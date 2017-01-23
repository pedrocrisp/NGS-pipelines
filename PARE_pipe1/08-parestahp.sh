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

#Args
sample=$1
bam_folder=$2
gff_ref=$3
genelist=$4
outfolder=$5

#Location of sample bam folder.
sample_dir="$bam_folder/$sample"

#Locate sample bam
sample_bam="$sample_dir/${sample}.bam"

#Create sample output directory. 
outdir="$outfolder/${sample}"
mkdir ${outdir}
results_table="$outdir/${sample}_parestahp.tab"

#Run parestahp
#USAGE: parestahp [options] -g GFF_FILE BAMFILE
#OPTIONS:
#    -u INT      Count INT bases upstream of the stop [default: 100]
#    -d INT      Count INT bases downstream stream of the stop [default: 100]
#    -g GFFFILE  GFF file describing gene models.
#    -l GENIDS   File containing a list of gene models. If not given, all gene
#                models are used, which can create inaccuate results. Please
#                provide a list of representative gene models.

parestahp \
-g $gff_ref \
-l $genelist \
$sample_bam \
> $results_table
