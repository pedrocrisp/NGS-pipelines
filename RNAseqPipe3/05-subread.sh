#!/bin/bash
set -e
set -x


###
#code to make it work on osx and linux
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

###
#reference sequence directory variable - user should create a link called subread_refdir in script dir to the location of the directory containing the subread indexfiles which must have the prefix "TAIR10_gen_chrc" (chrc means we included all 7 chromosomes)
refdir=$scriptdir/subread_refdir
#

sample=$1
sample_dir=reads_scythe_seqtk/$sample
outdir="align/${sample}"
mkdir ${outdir}
 
fastqs="$(ls $sample_dir/*.fq)"
numFqFiles=$(echo $fastqs | wc -w)

outsam="${outdir}/${sample}.sam"
outbam="${outdir}/${sample}" # no .bam, as samtools sort -f has a bug.
tmpbam="${outdir}/${RANDOM}.bam"

if [ ${numFqFiles} -eq 1 ]
then
echo subread-align -i ${refdir}/TAIR10_gen_chrc -r $fastqs -o "$outsam"
subread-align -i ${refdir}/TAIR10_gen_chrc -r $fastqs -o "$outsam"
elif [ ${numFqFiles} -eq 2 ]
then
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
echo subread-align -i ${refdir}/TAIR10_gen_chrc -r ${fq1} -R ${fq2} -o "$outsam"
subread-align -i ${refdir}/TAIR10_gen_chrc -r ${fq1} -R ${fq2} -o "$outsam"
else
echo "ERROR: not able to align multiple fq files per pair"
echo "fastqs:"
echo "${fastqs}"
exit 1
fi

echo "samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}"

samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}
