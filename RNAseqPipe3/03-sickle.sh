#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced. 
set -x

#Defines the sample that we are working with to the command line as the first token.
sample=$1
encoding=$2
read_ends=$3

#Specifies the directory in which the samples are located. 
inputDir=reads_noadapt

#Specifies the directory that the sample will be opened from.
sample_dir="${inputDir}/$sample"

#Defines the name of the output directory as 'reads_noadapt_trimmed'
outDir=reads_noadapt_trimmed

# Condition statement: if library is paired-end $3 == PE and trimming is done for both files.  If library is single end $3 == SE and single file is analysed.  If no  info is specified, script will error
if [ "$read_ends" == "PE" ]
then

#List all files ending with 'R1.noadapt.fq.gz' that are located within the specified sample directory and save these as the variable 'forward_fq.'
forward_fq="$(ls $sample_dir/*R1*.noadapt.fq.gz)"

#List all files ending with R2.noadapt.fq.gz that are located within the specified sample directory and save these as the variable 'reverse_fq.'
reverse_fq="$(ls $sample_dir/*R2*.noadapt.fq.gz)"

#Create the output directory and within this folder create a directory for the given sample. 
mkdir "${outDir}/$sample"

#Define the forward_fqname to be the name of the basic name of the forward_fq file (so get the name of the file rather than its full path). 
forward_fqname="$(basename $forward_fq)"

#Define forward_fq_outputfile path to be "reads_noadapt_trimmed/sample/forward_fqname.trimmed.fq'
forward_fq_outputFile="${outDir}/$sample/${forward_fqname%%.*}.trimmed.fq.gz"

#Define the reverse_fqname to be the name of the basic name of the reverse_fq file (so get the name of the file rather than its full path). 
reverse_fqname="$(basename $reverse_fq)"

#Define reverse_fq_outputfile path to be "reads_noadapt_trimmed/sample/reverse_fqname.trimmed.fq'
reverse_fq_outputFile="${outDir}/$sample/${reverse_fqname%%.*}.trimmed.fq.gz"

#Run the sickle program. 
#-f: specify forward read file for sample
#-r: specify reverse read file for sample
#-t: type of sequencing instrument
#-o: specifies location of output file for forward reads.
#-p: specifies location of output file for reverse reads. 
#-s: specifies location of output file for unpaired reads.
#-q: Set quality control. An averagq quality below 20 will be trimmed.  
#-l: Set length threshold in bp. Lengths below threshold are discarded. 
#-g for .gz input/output i suppose?
sickle pe \
-g \
-f $forward_fq \
-r $reverse_fq \
-t $encoding \
-o $forward_fq_outputFile \
-p $reverse_fq_outputFile \
-s "${outDir}/$sample/${sample}.trimmed.singles.fq.gz" \
-q 20 \
-l 20

elif [ "$read_ends" == "SE" ]
then

#List all files ending with 'R1.noadapt.fq.gz' that are located within the specified sample directory and save these as the variable 'forward_fq.'
forward_fq="$(ls $sample_dir/*.noadapt.fq.gz)"

#Create the output directory and within this folder create a directory for the given sample. 
mkdir "${outDir}/$sample"

#Define the forward_fqname to be the name of the basic name of the forward_fq file (so get the name of the file rather than its full path). 
forward_fqname="$(basename $forward_fq)"

#Define forward_fq_outputfile path to be "reads_noadapt_trimmed/sample/forward_fqname.trimmed.fq'
forward_fq_outputFile="${outDir}/$sample/${forward_fqname%%.*}.trimmed.fq.gz"

#Run the sickle program. 
#-f: specify forward read file for sample
#-r: specify reverse read file for sample
#-t: type of sequencing instrument
#-o: specifies location of output file for forward reads.
#-p: specifies location of output file for reverse reads. 
#-s: specifies location of output file for unpaired reads.
#-q: Set quality control. An averagq quality below 20 will be trimmed.  
#-l: Set length threshold in bp. Lengths below threshold are discarded. 
#-g for .gz input/output i suppose?
sickle se \
-g \
-f $forward_fq \
-t $encoding \
-o $forward_fq_outputFile \
-s "${outDir}/$sample/${sample}.trimmed.singles.fq.gz" \
-q 20 \
-l 20

else
echo "ERROR: SE or PE has not been specificed"

exit 1
fi
