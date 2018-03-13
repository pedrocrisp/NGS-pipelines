#!/bin/bash
#PBS -N 05-subread
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
module load samtools/0.1.18

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/05-runner.sh -j $threads -T $T -P $P -a $aligner -i $index -F $folder

######################
# to run eg
# 4 samples at once on 6 cores each (24 total)

# qsub -l walltime=2:00:00,nodes=1:ppn=24,mem=40gb \
# -v threads=4,T=6,P=3,aligner=subjunc,index=~/ws/refseqs/maize/subread_v1.6.0/Zea_mays.AGPv4,folder=reads_noadapt_trimmed \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/05-runner_qsub.sh

# bash /home/pete/gitrepos/NGS-pipelines/RNAseqPipe3/05-runner.sh-j 6 -T 1 -P 3 -a subjunc -i /home/pete/ws/refseqs/rice/rice_IR64_v7_MSU/subread/rice_ir64_v7
