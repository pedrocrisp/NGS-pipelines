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


