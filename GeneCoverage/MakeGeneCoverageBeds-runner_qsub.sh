#!/bin/bash
#PBS -N MakeGeneCoverageBeds
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
module load bedtools/2.25.0

bash ~/gitrepos/NGS-pipelines/GeneCoverage/MakeGeneCoverageBeds-runner.sh $alignFolder $threads $reference $coverage

######################
# to run eg
# 1 samples at once on 1 core

# qsub -l walltime=24:00:00,nodes=1:ppn=1,mem=50gb \
# -v alignFolder=tdf_for_igv,threads=1,reference=~/ws/refseqs/TAIR10/TAIR10_crisp_gene_primary.bed,coverage=full \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/GeneCoverage/MakeGeneCoverageBeds-runner_qsub.sh

# Make Full coverage beds for IGV viewing and my coverage plots
# bash ~/gitrepos/NGS-pipelines/GeneCoverage/MakeGeneCoverageBeds-runner.sh tdf_for_igv 1 ~/ws/refseqs/TAIR10/TAIR10_crisp_gene_primary.bed 5prime
