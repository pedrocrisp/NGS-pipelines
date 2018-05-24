#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash PerGeneCoverage_qsub.sh <sample_list.txt> <dataFolder> <filter_list_name> <gene_list_path> <library_layout>"

#define stepo in the pipeline - should be the same name as the script
step=GeneListMetaplotFunction

######### Setup ################
sample_list=$1
dataFolder=$2
filter_list_name=$3
gene_list_path=$4
library_layout=$5
if [ "$#" -lt "5" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for coverage analysis"
cat $sample_list
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t=1
else
qsub_t="1-${number_of_samples}"
fi
echo "argument to be passed to qsub -t is '$qsub_t'"

#find script to run, makes it file system agnostic
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi

########## Run #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make logs folder, timestamped
log_folder=logs/${timestamp}_${step}
mkdir $log_folder

#script path and cat a record of what was run
script_to_qsub=${scriptdir}/${step}.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub -t $qsub_t \
-o ${log_folder}/${step}_o \
-e ${log_folder}/${step}_e \
-v LIST=${sample_list},dataFolder=$dataFolder,filter_list_name=$filter_list_name,gene_list_path=$gene_list_path,library_layout=$library_layout \
$script_to_qsub

# to run
# bash /home/springer/pcrisp/gitrepos/NGS-pipelines/GeneCoverage2/GeneListMetaplotFunction_qsub.sh <sample_list.txt> <dataFolder> <filter_list_name> <gene_list_path> <library_layout>
# eg
# bash /home/springer/pcrisp/gitrepos/NGS-pipelines/GeneCoverage2/GeneListMetaplotFunction_qsub.sh samples.txt PerGeneCoverageBinned/per_gene_tables mRNAseq_RTLs ~/ws/refseqs/TAIR10/mRNAseq_RTLs.csv stranded
