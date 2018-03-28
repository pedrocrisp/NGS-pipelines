#!/bin/bash
#PBS -N 05b-samtools
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

######################
set -xeuo pipefail

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

mkdir -p logs

######################

module load parallel
# new modules
module load samtools/1.3_gcc-4.9.2_haswell

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/05b-runner.sh -j $threads -F $folder

######################
# to run eg
# 4 samples at once on 6 cores each (24 total)

# qsub -l walltime=12:00:00,nodes=1:ppn=24,mem=50gb \
# -v threads=4,folder=subjunc \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/05b-runner_qsub.sh

# sams to indexed and sorted bams from the subread-align folder
# bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/05b-runner.sh -j 4 -F subjunc
