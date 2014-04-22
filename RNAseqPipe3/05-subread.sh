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
#User defined reference sequence directory variable.
#This line defines the reference directory to be located in the script directory in a folder called 'subread_refdir.' The user should create a link for the subread_refdir in the script dir to map to the location of the directory containing their subread indexfiles. These index files can be created by using the 'subread-buildindex' program (refer to RNAseq pipeline user guide) and MUST have the prefix "TAIR10_gen_chrc" (chrc means we included all 7 chromosomes).
refdir=$scriptdir/subread_refdir
#

#Defines the sample that we are working with to the command line as the first token.
sample=$1

#Specifies the directory that the sample will be opened from.In this case, it is opening the sample folders located in the 'reads_noadapt_trimmed' folder.
sample_dir=reads_noadapt_trimmed/$sample

#Defines the output directory to be a folder with the sample name located within the 'align' directory. This will be used in the next step to create an output directory. 
outdir="align/${sample}"

#Creates an output directory to put the returned files to go into once subread has been run on the sample. The output from subread for the sample should go into a folder containing the sample's name, located within the 'align' directory.
mkdir ${outdir}

#List all files ending with 'trimmed.fq' that are located within the specified sample directory and save these as the variable 'fastqs.'
fastqs="$(ls $sample_dir/*trimmed.fq)"


numFqFiles=$(echo $fastqs | wc -w)

outsam="${outdir}/${sample}.sam"
outbam="${outdir}/${sample}" # no .bam, as samtools sort -f has a bug.
tmpbam="${outdir}/${RANDOM}.bam"

if [ ${numFqFiles} -eq 1 ]
then
echo subread-align -i ${refdir}/TAIR10_gen_chrc -r $fastqs -o "$outsam"
subread-align -i ${refdir}/TAIR10_gen_chrc -r $fastqs -o "$outsam"
elif [ ${numFqFiles} -eq 2 ]
then
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
echo subread-align -i ${refdir}/TAIR10_gen_chrc -r ${fq1} -R ${fq2} -o "$outsam"
subread-align -i ${refdir}/TAIR10_gen_chrc -r ${fq1} -R ${fq2} -o "$outsam"
else
echo "ERROR: not able to align multiple fq files per pair"
echo "fastqs:"
echo "${fastqs}"
exit 1
fi

echo "samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}"

samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}
