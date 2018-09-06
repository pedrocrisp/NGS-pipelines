#!/bin/bash
#PBS -N 05-bowtie2
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
module load bowtie2/2.3.4.1
module load samtools/1.3_gcc-4.9.2_haswell

bash ~/gitrepos/NGS-pipelines/PARE_pipe1/05-runner.sh $threads $reads $b_threads $reference $multimapping

######################
# to run eg 6 saamples in parallel with 3 cores each for mapping to maize B73 v4 genome

# qsub -l walltime=24:00:00,nodes=1:ppn=18,mem=80gb \
# -v threads=6,reads=reads_noadapt_cutadapt,b_threads=3,reference=~/ws/refseqs/maize/bowtie2/Zea_mays_AGPv4,multimapping=11 \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/PARE_pipe1/05-runner_qsub.sh
