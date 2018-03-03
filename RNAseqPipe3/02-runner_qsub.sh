#!/bin/bash -l
#PBS -l walltime=6:00:00,nodes=1:ppn=12,mem=40gb
#PBS -N 02-scythe
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

# Quality control on raw reads was performed with FASTQC v.0.11.5 (Andrews 2014)
# fastqc from RNAseq pipeline
# default is running 12 cores

module load parallel

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner.sh $1 $2

# to run
# qsub ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner_qsub.sh <number of threads> <prior>
# qsub ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner_qsub.sh 12 0.01
