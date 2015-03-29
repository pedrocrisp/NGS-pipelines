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

#cutadapt
#note the program has changed since it publication
#https://cutadapt.readthedocs.org/en/stable/guide.html#cutadapt-s-output


#Defines the sample that we are working with to the command line as the first token.
sample=$1
min=$2
max=$3
error_rate=$4

#Specifies the directory that the sample will be opened from.
sample_dir=reads/$sample

#List all files ending with '.fastq.gz' that are located within the specified sample directory and save these as the variable 'fastqs.'
fastqs="$(ls $sample_dir/*.fastq.gz)"

#Creates a new directory called 'reads_noadapt' and within this a folder for the sample. This creates the directory to put the output from scythe into (next step).

outdir="reads_noadapt_cutadapt/$sample"
mkdir $outdir

## uncomment if trying to keep reads failing length filter (this caused probelms in the past)
outdir_discard="reads_noadapt_cutadapt/${sample}/discard"
mkdir $outdir_discard



#This command runs cutadapt. It says 'for the fastqs listed within the sample directory do the following:
#1) keep the file names but remove the file extensions
#2) Store the output files in the specified sample folder within the reads_noadapt directory. Store the file name as 'samplename.noadapt.fq.gz' 
#3) Run the cutadapt function
#4) The Bash shell will replace the $(<...) with the content of the given file. The config file contains the -a: The adaptor sequences
#5) -m minimum read length discard reads if trimmed to smaller than m
#6) -M max rad length to keep
#7) -O only trimm adapter is if 10 nt match (elliminates trimming a few bases due to random matches to adapter eg the first 3). This option doesnt matter too much as it has to be 10 any way or else the farg witll be too long and filtered out anyway.
# --too-short-output rather than removing reads that are too short, write them to the discard folder for inspection
# --too-long-output ditto above for too long
# --untrimmed-output ditto if no adapt found, this is only approprite for sequencing of short fragemnts using longer read technology ie 20nt PARE tags or sRNAs using 50 nt reads.  Thus we expect an adapter to be read on the end or else something went wrong.  Not appropriate for RNAseq!
# -o output file (reads_noadapt/given sample name)
# input file

# not in use:
# --discard-untrimmed \
# --too-short-output $outputFile_too_short \
# --too-long-output $outputFile_too_long \
#
# --info-file $infoFile \
# --cut $end_trim \

for fq in $fastqs
do
fqname="$(basename $fq)"
outputFile="$outdir/${fqname%%.*}.noadapt.fq.gz"
#outputFile_too_short="$outdir_discard/${fqname%%.*}.too_short.fq.gz"
#outputFile_too_long="$outdir_discard/${fqname%%.*}.too_long.fq.gz"
outputFile_untrimmed="$outdir_discard/${fqname%%.*}.no_adapt_found.fq.gz"
#infoFile="$outdir/${fqname%%.*}.noadapt.info.txt"

cutadapt \
$(<cutadapt.conf) \
-e $error_rate \
-m $min \
-M $max \
-O 10 \
--untrimmed-output $outputFile_untrimmed \
-o $outputFile \
$fq
done

# histogram of read lengths (cut adapt log file isnt always informative enough... double check too)
# adapted from http://onetipperday.blogspot.com.au/2012/05/simple-way-to-get-reads-length.html
# textHistogram from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/textHistogram

echo "read length distribution for $outputFile"
zcat $outputFile | awk '{if(NR%4==2) print length($1)}' | textHistogram -maxBinCount=59 stdin
