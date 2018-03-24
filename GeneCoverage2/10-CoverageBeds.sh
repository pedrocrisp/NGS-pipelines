#!/bin/bash -l
#PBS -l walltime=6:00:00,nodes=1:ppn=1,mem=50gb
#PBS -N MakeGeneCoverageBeds
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAYID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################

module load bedtools/2.25.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID
# set sample = ID (could update all "samples" to be "ID" but this hack is shorter...)
sample=$ID
########## Run #################

sample_dir=$alignFolder/$sample
outFolder="${outdir}/${sample}"
mkdir ${outFolder}

####### stranded module
if [ "$coverage" == "5prime" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk -F$"\\t" 'BEGIN { OFS = FS } { $2=$2 "\\t" $2 } 1' $sample_dir/$sample.plus5.bed > $outFolder/$sample.plus.real.bed
awk -F$"\\t" 'BEGIN { OFS = FS } { $2=$2 "\\t" $2 } 1' $sample_dir/$sample.minus5.bed > $outFolder/$sample.minus.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.plus.real.bed > $outFolder/$sample.plus_sorted.bed
sort -k1,1 -k2,2n $outFolder/$sample.minus.real.bed > $outFolder/$sample.minus_sorted.bed
#remove intermediate file
rm -rv $outFolder/*real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.plus_sorted.bed -b $reference > $outFolder/$sample.plus.dist.bed
closestBed -D "ref" -a $outFolder/$sample.minus_sorted.bed -b $reference > $outFolder/$sample.minus.dist.bed
#remove intermediate file
rm -rv $outFolder/*_sorted.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.plus.dist.bed > $outFolder/$sample.plus.dist.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.minus.dist.bed > $outFolder/$sample.minus.dist.1k.bed

elif [ "$coverage" == "full" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk -F$"\\t" 'BEGIN { OFS = FS } { $2=$2 "\\t" $2 } 1' $sample_dir/$sample.plus.bed > $outFolder/$sample.plus.real.bed
awk -F$"\\t" 'BEGIN { OFS = FS } { $2=$2 "\\t" $2 } 1' $sample_dir/$sample.minus.bed > $outFolder/$sample.minus.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.plus.real.bed > $outFolder/$sample.plus_sorted.bed
sort -k1,1 -k2,2n $outFolder/$sample.minus.real.bed > $outFolder/$sample.minus_sorted.bed
#remove intermediate file

###### re comment after
rm -rv $outFolder/*real.bed
######

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.plus_sorted.bed -b $reference > $outFolder/$sample.plus.dist.bed
closestBed -D "ref" -a $outFolder/$sample.minus_sorted.bed -b $reference > $outFolder/$sample.minus.dist.bed
#remove intermediate file

###### re comment after
rm -rv $outFolder/*_sorted.bed
######

#subset to 1kb
awk -F$"\\t" '$NF<1000 && $NF>-1000' $outFolder/$sample.plus.dist.bed > $outFolder/$sample.plus.dist.1k.bed
awk -F$"\\t" '$NF<1000 && $NF>-1000' $outFolder/$sample.minus.dist.bed > $outFolder/$sample.minus.dist.1k.bed

####### non-stranded module

elif [ "$coverage" == "nonstranded" ]
then
#make a real bed file by adding chromosome end co-ordinates (same as start because nt resolution)
awk -F$"\\t" 'BEGIN { OFS = FS } { $2=$2 "\\t" $2 } 1' $sample_dir/$sample.bed > $outFolder/$sample.real.bed

#sort
sort -k1,1 -k2,2n $outFolder/$sample.real.bed > $outFolder/$sample.sorted.bed

#remove intermediate file
rm -rv $outFolder/*real.bed

#get distance measures
#reference is TAIR10_crisp_gene_primary.bed - fancy exon bed file prepared by SRE
closestBed -D "ref" -a $outFolder/$sample.sorted.bed -b $reference > $outFolder/$sample.dist.bed

#remove intermediate file
rm -rv $outFolder/*sorted.bed

#subset to 1kb
awk -F$'\t' '$NF<1000 && $NF>-1000' $outFolder/$sample.dist.bed > $outFolder/$sample.dist.1k.bed

else
echo "ERROR: coverage type has not been specificed select 5prime or full"
exit 1
fi

# Then remove any files absolutely required from making the plot, Just keep the dist.1k.bed
rm -rv $outFolder/*dist.bed
