#!/bin/bash
#PBS -N 03-sickle
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
# sickle

module load parallel
module load cutadapt/1.8.1

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/03-runner.sh $threads $min $max $error_rate

######################
# to run eg

# qsub -l walltime=12:00:00,nodes=1:ppn=18,mem=40gb \
# -v threads=18,min=18,max=25,error_rate=0.1 \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/03-runner_qsub.sh
