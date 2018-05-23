#!/bin/bash
#PBS -N 07-bam_to_tdf_stranded
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
module load samtools/1.5
module load bedtools/2.25.0

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/07-runner.sh $folder $threads $strand $chromoSizes

######################
# to run eg
# 2 samples at once on 2 cores each (24 total)

# qsub -l walltime=24:00:00,nodes=1:ppn=2,mem=50gb \
# -v folder=subjunc,threads=2,strand=stranded_PE,chromoSizes=~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/07-runner_qsub.sh

# #Now make new nt resolution beds to make bigWigs for fast IGV viewing and for coverage analysis; stranded libraries. Bedtools v2.26.0; bedGraphToBigWig v 4; IGV Version 2.3.35
# bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/07-runner.sh subjunc 2 stranded_SE ~/ws/refseqs/rice/rice_IR64_v7_MSU/all_con_chromosome_sizes
