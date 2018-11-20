#!/bin/bash -l
#PBS -l walltime=8:00:00,nodes=1:ppn=1,mem=50gb
#PBS -N tiles_bigWigs
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAYID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################

module load R/3.3.2
module load samtools/1.7
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID


########## Run #################
# args
echo ID
echo alignFolder
echo multi_filter
echo chrc_sizes

######
# dev
ID=B73_S
alignFolder=tmp_align
multi_filter=11

ID_dir=$alignFolder/$ID
results_folder="split_align/${ID}"
mkdir -p ${results_folder}

mkdir -p tmp_out_dir
tmp_out_dir="tmp_out_dir/${ID}"
mkdir -p ${tmp_out_dir}
######

# split
# quality filter
samtools view -q 10 $ID_dir/$ID.bam > $tmp_out_dir/${ID}_q10.sam

#21 nt
awk 'length($10) == 21 || $1 ~ /^@/' $tmp_out_dir/${ID}_q10.sam > $tmp_out_dir/${ID}_q10_21.sam

#22 nt
awk 'length($10) == 22 || $1 ~ /^@/' $tmp_out_dir/${ID}_q10.sam > $tmp_out_dir/${ID}_q10_22.sam

#24 nt
awk 'length($10) == 24 || $1 ~ /^@/' $tmp_out_dir/${ID}_q10.sam > $tmp_out_dir/${ID}_q10_24.sam


######

# multimapping summary stats and filter
# this step creates an NH:i tag in R, summarises NH:i distribution and then removes reads with NH:i equal to multi_filter

R -f ~/gitrepos/NGS-pipelines/smallRNAPipe1/16-split_summarise.R \
      --args ${ID} ${tmp_out_dir} ${multi_filter} ${results_folder}

########## clean up #################
# remove tmp dirs
rm -rv $tmp_out_dir
