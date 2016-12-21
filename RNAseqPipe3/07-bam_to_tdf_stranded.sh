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
#reference sequence directory variable - user should create a link called TAIR10_gen_chrc.chrom.sizes in script dir that points to the TAIR10_gen_chrc.genome.sizes file, or put a copy in there (chrc means we included all 7 chromosomes).

chrc_sizes=$4

sample=$1
alignFolder=$2
strand=$3
sample_dir=$alignFolder/$sample
outdir="tdf_for_igv/${sample}"
mkdir ${outdir}

# Condition statement: if library is stranded $3 == stranded and bigWigs are made for each strand.  If library is nonstranded $3 == nonstranded and nonstrandspecific bigWig is made.  If no strand info is specified, script will error
if [ "$strand" == "stranded_PE" ]
then


#split F and R reads into plus and minus strand taking into account PE
#http://seqanswers.com/forums/showthread.php?t=29399

#R1 forward strand
samtools view -f 99 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1F.bam

#R2 reverse strand
samtools view -f 147 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2R.bam

samtools merge -f $sample_dir/${sample}.forward.bam $sample_dir/${sample}.R1F.bam $sample_dir/${sample}.R2R.bam

#R1 reverse strand
samtools view -f 83 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1R.bam

#R2 forward strand
samtools view -f 163 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2F.bam

samtools merge -f $sample_dir/${sample}.reverse.bam $sample_dir/${sample}.R1R.bam $sample_dir/${sample}.R2F.bam

rm $sample_dir/${sample}*.R*.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bed

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

bedGraphToBigWig $outdir/*.plus.bg $chrc_sizes $outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/*.minus.bg $chrc_sizes $outdir/$sample.minus.bigWig

elif [ "$strand" == "reverse_stranded_PE" ]
then


#split F and R reads into plus and minus strand taking into account PE
#http://seqanswers.com/forums/showthread.php?t=29399

#R1 forward strand
samtools view -f 99 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1F.bam

#R2 reverse strand
samtools view -f 147 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2R.bam

samtools merge -f $sample_dir/${sample}.forward.bam $sample_dir/${sample}.R1F.bam $sample_dir/${sample}.R2R.bam

#R1 reverse strand
samtools view -f 83 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R1R.bam

#R2 forward strand
samtools view -f 163 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.R2F.bam

samtools merge -f $sample_dir/${sample}.reverse.bam $sample_dir/${sample}.R1R.bam $sample_dir/${sample}.R2F.bam

rm $sample_dir/${sample}*.R*.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bed

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

bedGraphToBigWig $outdir/*.plus.bg $chrc_sizes $outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/*.minus.bg $chrc_sizes $outdir/$sample.minus.bigWig

elif [ "$strand" == "stranded_SE" ]
then

#R1 forward strand -f 0 does not work, all reads retained... "samtools view -f 0 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam"
#not 100% sure we want everything that is NOT unmapped or reserse mapped, not sure what else coudl be included... of well
samtools view -F 4 -b $sample_dir/$sample.bam | samtools view -F 16 -b - > $sample_dir/${sample}.forward.bam

#R1 reverse strand
samtools view -f 16 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.reverse.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bed

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

bedGraphToBigWig $outdir/*.plus.bg $chrc_sizes $outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/*.minus.bg $chrc_sizes $outdir/$sample.minus.bigWig

### 5' analysis
### Note consider adding this to other sections, alternatively leaving it here will probably mean it only gets run for PARE data (single end stranded)
#stranded bedgraphs coverage of 5' position only
bedtools genomecov -5 -d -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus5.bed
#minus strand reads bedgraph
bedtools genomecov -5 -d -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus5.bed

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -5 -bga -scale -1 -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.minus5.bg
#minus strand reads bedgraph
bedtools genomecov -5 -bga -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.plus5.bg

bedGraphToBigWig $outdir/*.plus5.bg $chrc_sizes $outdir/$sample.plus5.bigWig
bedGraphToBigWig $outdir/*.minus5.bg $chrc_sizes $outdir/$sample.minus5.bigWig
###########

##

elif [ "$strand" == "reverse_stranded_SE" ]
then

#R1 forward strand -f 0 does not work, all reads retained... "samtools view -f 0 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam"
#not 100% sure we want everything that is NOT unmapped or reserse mapped, not sure what else coudl be included... of well
samtools view -F 4 -b $sample_dir/$sample.bam | samtools view -F 16 -b - > $sample_dir/${sample}.forward.bam

#R1 reverse strand
samtools view -f 16 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.reverse.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -ibam $sample_dir/*reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam $sample_dir/*forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bed

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

bedGraphToBigWig $outdir/*.plus.bg $chrc_sizes $outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/*.minus.bg $chrc_sizes $outdir/$sample.minus.bigWig


##



elif [ "$strand" == "nonstranded" ]
then

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph


#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
#could add this option to all steps above - output .bed file is about the same size as the bam ie liek 1.8 GB... quite big!
bedtools genomecov -d -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bed

#make tdf
echo "bedgraph to binary tiled data (.tdf) file"
igvtools toTDF $outdir/*.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/*.bedgraph $chrc_sizes $outdir/$sample.bigWig

else
echo "ERROR: it has not been specificed whether library is stranded on not"

exit 1
fi


