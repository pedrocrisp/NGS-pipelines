#!/bin/bash
#PBS -N 04-fastqc
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
# fastqc post triming

module load parallel
module load fastqc/0.11.5

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/04-runner.sh $threads

######################
# to run eg

# qsub -l walltime=2:00:00,nodes=1:ppn=18,mem=40gb \
# -v threads=18 \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/04-runner_qsub.sh
