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
# scythe

module load parallel

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/03-runner.sh $encoding $threads $read_ends

######################
# to run eg

# qsub -l walltime=6:00:00,nodes=1:ppn=18,mem=40gb \
# -v encoding=sanger,threads=12,read_ends=PE \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/03-runner_qsub.sh
