#RNAseqPipe3 README

---

For detailed descriptions of the tools implemented in this pipeline see [documentation by Bec](https://github.com/pedrocrisp/NGS-pipelines/blob/master/Docs/RNAseq.md)

---

Pipeline to process RNAseq data

Reads must be in a folder called reads, with subdirectories for each sample.  Each sample directory should contain one fasatq.gz; if there are multiple fastq files cat them together using 10-runner.sh (10-cat.sh) called as below.

Scripts are called from the parent directory of the reads folder

Each script and std error is printed to the logs directory

Each step is called in order and run in bash eg 

```
bash ~/path_to_01.runner.sh
```

To run the whole pathway run desired scripts in order, eg the following 

1. Cats together the 2 files that come off the seqeunce for each sample
2. Runs fastQC
3. Trims adapters
4. Quality trims
5. Runs fastQC again
6. Aligns 3 samples in parallel with subread using 2 cores per sample
7. Compresses the sam file output of subread into sorted and indexed bam for use in IGV etc on 9 cores in the folder subread-align
8. Summarise reads to features using featureCounts, set strandedness to 2 for stranded dUTP method (2 for Illumina stranded kit), and find the bams in the folder subread-align
8. Makes tdf files and strand seperated bigWig files from the bams for quicker IGV loading using bams in subread-align folder

```
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/10-runner.sh;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/01-runner.sh;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/02-runner.sh;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/03-runner.sh;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/04-runner.sh;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/05-runner.sh -j 3 -P 2 -a subread-align;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/05b-runner.sh -j 9 -F subread-align;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/06-runner.sh 2 subread-align;
bash ~/gitrepos/ngs-pipelines/RNAseqPipe3/07-runner.sh subread-align
```

---
##Step 01-fastqc.sh
```
usage:  01-runner.sh <number of threads passed to parallel-j flag>
```

---
##Step 02-scythe.sh - adapter removal

```
usage:  02-runner.sh <number of threads passed to parallel-j flag>
```

Step 02-scythe.sh requires the file "truseq_adapters.fasta" to be in the script directory or a symbolic link to it

---
##Step 03-sickle.sh -quality trimming
```
usage:  03-runner.sh <number of threads passed to parallel-j flag>
```
Trims until average quality is phred 20 in sliding window

---
##Step 04-fastqc.sh
```
usage:  04-runner.sh <number of threads passed to parallel-j flag>
```

---
##Step 05-subread.sh -genome alignment
```
usage:  05-runner.sh [-j value] [-P value] [-a <aligner>]
```
Must provide 
1. -j the number of jobs/samples to run in parallel 
2. -P the number of threads subread-align/subjunc uses per sample
3. -a whether to use subread-align or subjunc 

Step 05-subread.sh also requires a folder (or symbolic link) called "subread\_refdir" to be in the same folder as the script. The "subread\_refdir" folder contains the index files, those files must have the prefix "TAIR10\_gen\_chrc". Make sure when you build the index that the headers of the chloroplast and mitchondria fasta files agree with other reference files eg the chromosome sizes file.  I channge the headers to ChrC and ChrM (note the captials), this oversight causes days of frustration...

```
mkdir subread.v1.4.5
cd subread.v1.4.5
subread-buildindex -o TAIR10_gen_chrc ../chromosomes/TAIR10_chr*
#add the .saf file to this folder
cp ...TAIR10_GFF3_genes.saf ./
```

*NOTE: percy inexplicably cant handle more than 1 (-j 1) parallel job at once at this step, memory errors evrytime on one or 2 samples

---
##Step 05b-samtools.sh - sam to indexed, sorted bam
```
usage:  05b-samtools.sh [-j value] [-a <aligner>]
```
Must provide 

1. -j the number of jobs/samples to run in parallel 
2. -a whether to use subread-align or subjunc 

---
##Step 06-featureCounts.sh - read sumarisation
```
usage:  06-runner.sh <strandedness> <alignment_folder> <threads>
```

Step 06-featureCounts.sh requires the "subread\_refdir" to also contain a TAIR10\_gen_chrc.chrom.sizes file

---
##Step 07-bam\_to\_tdf_stranded.sh - make bigWigs etc

```
usage:  07-bam_to_tdf_stranded.sh <alignment_folder> <threads> <strandedness of library>
```

Strandedness can be "stranded_PE", "stranded_SE" or "nonstranded".

---
Link ref seq on the server:

```
ln -s /home/pete/workspace/refseqs/TAIR10_gen/ subread_refdir

ln -s /home/pete/workspace/refseqs/TAIR10_gen/TAIR10_gen_chrc.chrom.sizes TAIR10_gen_chrc.chrom.sizes

ln -s /home/pete/workspace/Illumina/truseq_adapters.fasta truseq_adapters.fasta

```

## Example

```
############
#Quality control on raw reads was performed with FASTQC v.0.11.2 (Andrews 2014)

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/01-runner.sh 20

############
#Adapters were removed using scythe v.0.991 with flags -p 0.01 for the prior (Buffalo 2014)
#Contamination rate ~0.03%

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/02-runner.sh 12

############
#Reads were quality trimmed with sickle v.1.33 with flags -q 20 (quality threshold) -l 20 (minimum read length after trimming) (Najoshi 2014)

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/03-runner.sh 20

############
#Quality control on filtered reads was performed with FASTQC v.0.11.2 (Andrews 2014)

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/04-runner.sh 20


############
#align, use subjunc (if on Percy, one one job at a time or it dies...)
#The trimmed and quality filtered reads were aligned to the Arabidopsis genome (TAIR10) using the subjunc v.1.4.6 aligner with -u and -H flags to report only reads with a single unambiguous best mapping location (Liao, Smyth, and Shi 2013).

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/05-runner.sh -j 12 -P 1 -a subjunc

############
# Reads were then sorted, indexed and compressed using samtools v1.1-26-g29b0367 (\cite{li_sequence_2009}) 

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/05b-runner.sh -j 9 -F subjunc

############
#featureCounts summarise counts reverse stranded (dUTP)
#For standard differential gene expression testing the number of reads mapping per TAIR10 gene loci was summarised using featureCounts v. 1.4.6 with flags -p and -C to discard read pairs mapping to different chromosomes and the -s flag set to 0 for a non-strand specific library, multimapping reads and multioverlapping reads (reads mapping to overlapping regions of more than one gene loci) were not counted (Liao, Smyth, and Shi 2014). Reads were summarised to parent TAIR10 gene loci rather than individual splice variants by summarising to the genomic coordinates defined by the feature "gene" in the TAIR10_GFF3_genes.gff reference (last modified 14/10/2010 ftp://ftp.arabidopsis. org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff).

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/06-runner.sh 0 subjunc 12 TAIR10_GFF3_genes.saf

############
#make bigWigs for fast IGV viewing
# Strand specific bigwig files were generated using bedtools genomecov v2.16.1 (\cite{quinlan_bedtools:_2010}) and the UCSC utility bedGraphToBigWig \href{http://hgdownload.cse.ucsc.edu/admin/exe/}{\textit{link}} for viewing in IGV (\cite{robinson_integrative_2011}).

bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/07-runner.sh subjunc 12 nonstranded

############


```
