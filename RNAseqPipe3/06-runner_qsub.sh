#!/bin/bash
#PBS -N 06-featureCounts
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

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/06-runner.sh $strand $alignFolder $threads $reference

######################
# to run eg
# reverse stranded (Illumina strand specific), subjunc folder, 12 threads, SAF reference?

# qsub -l walltime=6:00:00,nodes=1:ppn=24,mem=40gb \
# -v strand=2,alignFolder=subjunc,threads=12,reference=~/ws/refseqs/maize/Zea_mays_AGPv4_36_gene.SAF \
# -o logs \
# -e logs \
#  ~/gitrepos/NGS-pipelines/RNAseqPipe3/06-runner_qsub.sh

#featureCounts summarise counts Illumina stand-specific (reverse stranded library)
#For standard differential gene expression testing the number of reads mapping per IR64_v7 gene loci was summarised using featureCounts v. 1.5.0-p1 with flags -p and -C to discard read pairs mapping to different chromosomes and the -s flag set to 2 for a strand specific library, multimapping reads and multioverlapping reads (reads mapping to overlapping regions of more than one gene loci) were not counted (Liao, Smyth, and Shi 2014). Reads were summarised to parent IR64_v7 gene loci rather than individual splice variants by summarising to the genomic coordinates defined by the feature "gene" in the all.gff reference (last modified 7/2/2012 ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3).

#80-87% of mapped reads mapped to features. alx8 samples were 79.0,81.6 and 80.2, xrn 83.6, 84.2 and 84.9, WT 87.3, 87.9 and 86.6. SO it appears that the mutants have more non-coding RNA (assuming that this isnt just a systematic increase in rRNA contamination, which is a possibility that would have to be addressed later).

# bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/06-runner.sh 2 subjunc 9 ~/ws/refseqs/rice/rice_IR64_v7_MSU/featureCounts/all.saf
