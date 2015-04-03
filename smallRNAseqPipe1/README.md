#sRNAseqPipe1 README

---

For detailed descriptions of some of the tools implemented in this pipeline see [documentation by Bec](https://github.com/pedrocrisp/NGS-pipelines/blob/master/Docs/RNAseq.md)

---

Pipeline to process sRNAseq data

Reads must be in a folder called reads, with subdirectories for each sample 

Scripts are called from the parent directory of the reads folder

Each script and std error is printed to the logs dircetory

Each step is called in order and run in bash eg 

See scripts for specific usage and required args, all are sent to parallel with an arg to set the number of parallel job (threads).

```
bash ~/path_to_01.runner.sh
```

To run all the main pathway script call the 0*-runner.sh at once, eg:

```
for script in ~/<path_to_pipeline_folder>/0*runner.sh ; do echo $script ; bash $script ; done
```

---

Steps 02-scythe.sh has worked with much sucess (due ot the high level of adapter contamination adn low read quality at the 3' end?) requires the file "truseq_adapters.fasta" to be in the script directory or a symbolic link to it

---

Steps 03-cutadapt.sh requires the file "cutadapt.conf" which contains the adapter option and sequence (eg -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC) to be in the working directory or a symbolic link to it

For NEB Next libraries -a AGATCGGAAGAGC

It is recommended to filter reads by length -m 18, -M 25 after trimming.  Shorter reads take forever to align if report all alignments is selected for bowtie2.

---

Steps 05-bowtie2.sh requires a symbolic link called "bowtie2_refdir" to the folder containing the index files, those files must have the prefix "TAIR10_allchr"

On the server:

```
ln -s /home/pete/workspace/refseqs/TAIR10_gen/ bowtie2_refdir

ln -s /home/pete/workspace/Illumina/truseq_sRNA_adapters.fasta truseq_adapters.fasta

```

## Example

```
########### pipeline.sh ##########
#This is the ultimate analysis pipeline for RRGS sRNA data. All scripts called from the parent directory of the project folder, which contains a sub folder called reads, containing a folder for each sample, which in trun contains the fastq files.

###
#Fastq files were provided as 2 files by the sequencing facility, denoted by an L001 and L002, these were concatinated together and renamed.

for FILE in *_317* ; do mv -v "$FILE" "${FILE//Raw/Sample}" ; done

bash ~/gitrepos/NGS-pipelines/tools/cat_fastqs-runner.sh 11 reads SE

###
#Read quality control was performed with FastQC v0.11.2 (\cite{andrews_fastqc_????}) before and after trimming.

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/01-runner.sh 11 reads

###
# Given that reads were sequenced with 50 bp reads length, all sRNAs were expected to contain adpater on the 3'. Adapters were trimmed using cutadapt v1.8 (Python2.7.8; \cite{martin_cutadapt_2011}; note that the cutadapt algorithm has changed significantly since publication \href{https://cutadapt.readthedocs.org/en/stable}{link}). Cutadapt was run with the flags: -e 0.1, first 9 bases of adapter must match perfectly; -m 18 -M 25, only keep reads lengths from 18 nt to 25 nt after timming; -O 10, adapter must overlap by at least 10 bases; -a AGATCGGAAGAGC, adapter sequence. Histograms of read-length distribution were calculated using textHistogram \href{http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/textHistogram}{\textit{link}}.

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/03-runner.sh 20 18 25 0.1

###
# fastqc

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/01-runner.sh 11 reads_noadapt_cutadapt_18_25

###
# Reads were then aligned to the Arabidopsis genome (TAIR10) using bowtie2 v2.2.5 (\cite{langmead_fast_2012}), using the flags: -a to report all matches for multi-mapped reads; -D 20 and -R 3, increases the likelihood that bowtie2 will report the correct alignment for a read that aligns many places and -i S,1,0.50, reduces the substring interval, further increasing sensitivity; --end-to-end, preventing trimming of reads to enable alignment; -L 10 reduces substring length to 10 (default 22) as these are short reads and -N 0 requires exact match in the seed; --score-min L,0,0 reports only exact matches in --end-to-end mode (alignment score of 0 required which is max possible in end mode).

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/05-runner.sh 2 reads_noadapt_cutadapt 6

###
# Reads were then sorted, indexed and compressed using samtools v1.1-26-g29b0367 (\cite{li_sequence_2009}) and strand specific bigwig files were generated using bedtools genomecov v2.16.1 (\cite{quinlan_bedtools:_2010}) and the UCSC utility bedGraphToBigWig \href{http://hgdownload.cse.ucsc.edu/admin/exe/}{\textit{link}} for viewing in IGV (\cite{robinson_integrative_2011}).
bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/07-runner.sh align_bowtie2 12 stranded_SE

###########################################




```
