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
#This line specifies that the reference directory is located in the script directory in a folder called 'subread_refdir.' The user should create a link for the subread_refdir in the script dir to map to the location of the directory containing their subread saf file. The subread saf file must be called "TAIR10_GFF3_genes.saf"
refdir=$scriptdir/subread_refdir
#

#Defines the sample that we are working with to the command line as the first token.
sample=$1
alignFolder=$3
SAF=$4

#Specifies the directory that the sample will be opened from. In this case, it is opening a sample folder located in the 'align' folder.
sample_dir=$alignFolder/$sample

#Defines the output directory to be a folder with the sample name located within the 'featureCounts' directory. This will be used in the next step to create an output directory.
outdir="featureCounts/${sample}"

#Creates an output directory to put the returned files to go into once featureCounts has been run on the sample. In this case, the output from featureCounts for the sample should go into a folder containing the sample's name, located within the 'featureCounts' directory.
mkdir ${outdir}

############### potentially counting multimapping reads
# featureCounts on all reads in bam - not this will count multimapping reads from bowtie multiple times

#This command runs featureCounts.
#-F: Specify the format for the annotated genome file you are using.
# ##NOT USED## -p: For paired end reads only, specified fragments will be counted as opposed to reads.
# ##NOT USED## -C: Specifies that chimeric fragments (those fragments with their two ends aligned to different chromosomes) will NOT be counted.(Used in conjunction with -p).
#-s: Specifies if strand specific read counting should be performed. Three values are possible: o =unstranded, 1= stranded, 2= reversely stranded.
#-a: specify the filepath to your annotated genome library.
#-o: Output directory.
#${sample_dire}/*.bam - specifies that the input for feature counts is located in the named sample directory.

strand=$2

featureCounts\
    -F SAF\
    -s $strand\
    -a $SAF\
    -o "$outdir/${sample}.counts"\
    "${sample_dir}/${sample}.bam"


    ############### filter to unique reads to omit counting putative multimapping reads
    # note this omits all reads with secondary alignments, some may not be as good as the primary but the primary is omited anyway
    # featureCounts

    #This command runs featureCounts.
    #-F: Specify the format for the annotated genome file you are using.
    # ##NOT USED## -p: For paired end reads only, specified fragments will be counted as opposed to reads.
    # ##NOT USED## -C: Specifies that chimeric fragments (those fragments with their two ends aligned to different chromosomes) will NOT be counted.(Used in conjunction with -p).
    #-s: Specifies if strand specific read counting should be performed. Three values are possible: o =unstranded, 1= stranded, 2= reversely stranded.
    #-a: specify the filepath to your annotated genome library.
    #-o: Output directory.
    #${sample_dire}/*.bam - specifies that the input for feature counts is located in the named sample directory.

    # remove reads with secondary alignments
    # samtools view Sample_WT_277_1.bam -h | grep -E "@|NM:" | grep -v "XS:" | samtools view -h -b -o Sample_WT_277_1_uniq.bam -
    samtools view -h "${sample_dir}/${sample}.bam" | grep -E "@|NM:" | grep -v "XS:" | samtools view -h -b -o "${sample_dir}/${sample}_uniq.bam" -

    featureCounts\
        -F SAF\
        -s $strand\
        -a $SAF\
        -o "$outdir/${sample}_uniq.counts"\
        "${sample_dir}/${sample}_uniq.bam"
