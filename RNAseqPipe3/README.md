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
usage:  06-runner.sh [-j value] [-P value] [-a <aligner>]
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
usage:  06-runner.sh <strandedness> <alignment_folder>
```

Step 06-featureCounts.sh requires the "subread\_refdir" to also contain a TAIR10\_gen_chrc.chrom.sizes file

---
##Step 07-bam\_to\_tdf_stranded.sh - make bigWigs etc
```
usage:  07-bam_to_tdf_stranded.sh <alignment_folder>
```
---
Link ref seq on the server:

```
ln -s /home/pete/workspace/refseqs/TAIR10_gen/ subread_refdir

ln -s /home/pete/workspace/refseqs/TAIR10_gen/TAIR10_gen_chrc.chrom.sizes TAIR10_gen_chrc.chrom.sizes

ln -s /home/pete/workspace/Illumina/truseq_adapters.fasta truseq_adapters.fasta

```

