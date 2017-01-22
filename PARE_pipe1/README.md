# PARE analysis Pipeline

This Pipeline processes PARE (degradome) data from fastq, adapter removal to genome mapping and summarisation. It does not perform sRNA target prediction, rather is designed for interrogation of the RNA degradome beyond miRNA targets.

## QC -fastqc

FastQC is used to examine the fastq files.  Typically sequencing runs into the 3' adapter because in this protocol 20-12 nt fragemnts are cloned yet seqeuncing is usuallt longer (50-100 bp).  Hence, data quality is often terrible beyond base 21 due to colour imbalance on the sequencer.

This script can now be run on fastqs before and after trimming by specifying the reads folder. Also it recognises fasq.gz and fq.gz (searches for *.q.gz).

usage="USAGE:

```01-runner.sh <number of threads> <reads folder>"```

Number of threads refers to parallel ie how many samples to run at once (the rest are queued)


## Adapter removal -cutadapt

Adapter removal can present a challenge because of the colour imbalance owing to the extremely consistent 20-21 nt fragemnt size. The approach I have settled on for the moment is to simply cut the reads back to 20 nt (21 nt reads lose a base but seemed better than including an adapter base on the 20 nt reads, this way can map 20 nt reads accepting perfect matches only). The 20 nt reads are checked for adapters, (3 nt perfect match min) but very few reads are filtered out at this stage. Using this appraoch >90% of reads should map (80% uniquely), so seems satisfactory. The distribution of read lengths is also calculated folowwing cutadapt and written into the log file, requires [textHistogram](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/textHistogram), module adapted from [One tip per day](http://onetipperday.blogspot.com.au/2012/05/simple-way-to-get-reads-length.html).

usage="USAGE:

```03-runner.sh <number of threads> <min length> <max length> <error rate> <3 end trim length>"```

Positional args

1. Number of threads refers to parallel ie how many samples to run at once (the rest are queued)
2. min fragment length to retain after trimming
3. max fragment length to retain after trimming
4. 4. error rate - default is 0.1 (10%)*match length, recommended sticking with this ie 1 error per 10 base match
5. bases to trim from end, negative number indicates trim from 3' end, positive number trim 5' - recommended to trim back to 20 nt (-31 for 50 bp run, -81 for 100 bp run)

cutadapt requires a file called cutadapt.conf to be in the directory that the script is called from.  The file is subsitituted as an argument to cutadapt as contains the adapter info.  For these PARE libraries the adapter is the 3' TruSeq sRNA adpater:

```
# > PARE Final Product
# AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACNNNNNNNNNNNNNNNNNNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
# TTACTATGCCGCTGGTGGCTGTCCAAGTCTCAAGATGTCAGGCTGNNNNNNNNNNNNNNNNNNNNACCTTAAGAGCCCACGGTTCCTTGAGGTCAGTGNNNNNNTAGAGCATACGGCAGAAGACGAAC
```

So this file should read ``` -a <adapter sequence>```:

```
cat cutadapt.conf 
-a TGGAATTCTC
```

Note: there is a scythe script in here too but scythe kept trimming almost everything away leaving no reads, couldnt git it too work, maybe it is confused by  the low quality of reads, cutadapt worked well enough so didnt persist, although scythe may well do a better job if it could be optimized.

## Alignment - bowtie2

For alignment bowtie2 is used to map to the TAIR10 genome, requiring perfect match (because these are short reads...). 

`05-runner.sh <number of threads> <reads folder> <bowtie threads per job> <reference>`

## Example:

```
########### pipeline.sh ##########i
#This is the ultimate analysis pipeline for RRGS PARE data. All scripts called from the parent directory of the project folder, which contains a sub folder called reads, containing a folder for each sample, which in trun contains the fastq files.

############ #Step 1 cat files ############
#Step 1 cat files
#Fastq files were provided as 2 files by the sequencing facility, denoted by an L001 and L002, these were concatinated together and renamed.
bash ~/gitrepos/NGS-pipelines/tools/cat_fastqs-runner.sh 11 reads SE

########### #Step 2 QC raw fastq ###########
#Step 2 QC raw fastq
#Quality control on raw reads was performed with FASTQC v.0.11.2 (Andrews 2014)
#fastqc from RNAseq pipeline
bash ~/gitrepos/NGS-pipelines/RNAseqPipe3/01-runner.sh 20
#note, since i have initiated a PARE pipeline with a new fastqc script (specify read dir).

########### #Step 3 trim adapters ###########
# Step 3 trim adapters
# Reads were trimmed using cutadapt v1.8 (Python2.7.8; \cite{martin_cutadapt_2011} however the cutadapt algorithm has change significantly since publication \href{https://cutadapt.readthedocs.org/en/stable
# Adapters are TruSeq sRNA adapters, 3' end as per sRNA protocol, 5' slightly shorter adapter (hence the custom sequencing primer)
# > PARE Final Product
# AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACNNNNNNNNNNNNNNNNNNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
# TTACTATGCCGCTGGTGGCTGTCCAAGTCTCAAGATGTCAGGCTGNNNNNNNNNNNNNNNNNNNNACCTTAAGAGCCCACGGTTCCTTGAGGTCAGTGNNNNNNTAGAGCATACGGCAGAAGACGAAC
# Strategy:
# cut reads back to 20 nt
# search for TGGAATTCTC using cutadapt
# filter short reads too (keep reads 14-22 nt)

bash ~/gitrepos/NGS-pipelines/PARE_pipe1/03-runner.sh 11 14 22 0.1 -31

#fastqc again
bash ~/gitrepos/NGS-pipelines/PARE_pipe1/01-runner.sh 11 reads_noadapt_cutadapt

########### #Step 4 align to genome ###########
# Reads were then aligned to the Arabidopsis genome (TAIR10) using bowtie2 v2.2.5 (\cite{langmead_fast_2012}), using the flags: -a to report all matches for multi-mapped reads; -D 20 and -R 3, increases the likelihood that bowtie2 will report the correct alignment for a read that aligns many places and -i S,1,0.50, reduces the substring interval, further increasing sensitivity; --end-to-end, preventing trimming of reads to enable alignment; -L 10 reduces substring length to 10 (default 22) as these are short reads and -N 0 requires exact match in the seed; --score-min L,0,0 reports only exact matches in --end-to-end mode (alignment score of 0 required which is max possible in end mode).

bash ~/gitrepos/NGS-pipelines/PARE_pipe1/05-runner.sh 11 reads_noadapt_cutadapt_20nt 2

########### #Step 4 bam to wig files ###########
# Reads were then sorted, indexed and compressed using samtools v1.1-26-g29b0367 (\cite{li_sequence_2009}) and strand specific bigwig files were generated using bedtools genomecov v2.16.1 (\cite{quinlan_bedtools:_2010}) and the UCSC utility bedGraphToBigWig \href{http://hgdownload.cse.ucsc.edu/admin/exe/}{\textit{link}} for viewing in IGV (\cite{robinson_integrative_2011}).

bash ~/gitrepos/NGS-pipelines/smallRNAseqPipe1/07-runner.sh align_bowtie2 12 stranded_SE

```

