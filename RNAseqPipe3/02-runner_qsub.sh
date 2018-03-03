#!/bin/bash -l
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
# scythe

module load parallel

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner.sh $threads $prior

######################
# to run eg

# qsub -l walltime=6:00:00,nodes=1:ppn=18,mem=40gb \
# -v threads=18,prior=0.01 \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner_qsub.sh
