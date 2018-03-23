#!/bin/bash
#PBS -N GeneCoverageR
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
module load R/3.3.2

bash ~/gitrepos/NGS-pipelines/GeneCoverage/GeneCoverageR-runner.sh $alignFolder $threads $library_layout

######################
# to run eg
# 1 samples at once on 1 core

# qsub -l walltime=24:00:00,nodes=1:ppn=1,mem=50gb \
# -v alignFolder=tdf_for_igv_coverage_beds,threads=1,library_layout=nonstranded \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/GeneCoverage/GeneCoverageR-runner_qsub.sh

#To run:
#bash ~/gitrepos/NGS-pipelines/GeneCoverage/GeneCoverageR-runner.sh tdf_for_igv_coverage_beds 1
