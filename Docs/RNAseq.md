RNAseq pipeline user guide
===========================

Fast-QC
-------

**Function:**

FastQC is a java based program designed to perform a quality control analysis on your reads. This is an important tool as it enables you to identify any issues that may have occurred during high throughput sequencing (eg. decrease in read quality towards the end of every read) or any issues with your starting library material.  

**Use:** 

Fast QC can be used to perform simple quality control analysis on your raw data, and can also be used after various steps in this pipeline (eg. removal of adapter sequences) to see if these programs have improved read quality. 

**Requirements:**

Fast QC can read files ending in: 
* FastQ (all quality encoding variants), 
* Casava FastQ files^, Colorspace FastQ, 
* GZip compressed FastQ, 
* SAM, 
* BAM or 
* SAM/BAM Mapped only^ 

^Indicates that you will need to tell FastQC you are using SAM/BAM or Casava files, as FastQC will not automatically select these formats. 

**Evaluation of results:**

* Once you have run Fast-QC, the results report will be available as a http file. The report should open in your browser when you click on the file.
* You will see that the analysis is segmented into a series of modules, with a summary of the module outcomes provided on the left hand side. Ratings include: normal (green tick), slightly abnormal (orange exclamation mark) and very abnormal (red cross). 
* Here is an example of a [good quality](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc/fastqc_report.html) and [poor quality] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc/fastqc_report.html) read. A brief description of the purpose of each module included in a FastQC report is provided below. 

**Basic statistics module-** Provides a summary of the data, including; file name/type, total number of sequences processed (actual and estimated), range in sequence length (if all reads are the equal length, only one value will be reported) and overall GC content expressed as a percentage. 

*N.b: basic statistics will always appear to have a green tick (never flags a warning/ error), so further analysis of data should always occur*

**Per base sequence quality-** Shows the variation in sequence quality at every base position for all reads in the file. The y axis represents the quality scores and the x axis is the base position in the reads. The information is presented as a boxplot for each read position, with the red line representing the median and the blue line representing the mean values. As base calling quality increases with quality number, the graph has been divided into green, orange and red sections, respectively representing good, satisfactory and poor quality base calls. In a good quality run (refer to example shown above), one would want their boxplots to be in the green section and have little variation. 

**Per sequence quality scores-** Enables you to see if a subset of your reads have overall poor quality base recall. Often a small percentage of your reads will be universally poor quality, generally because they have been poorly imaged (eg. may be on edge of the field of view). If there is a significant percentage of reads that have poor quality in a region of base pairs, this may suggest a systematic problem. 

**Per base sequence content-** Provides a plot of the percentage of each base called (A, T, G, C) at each base pair position within the library file. In a random library, the proportion of the four bases for a given base pair position should theoretically be approximately equal, thus the lines for the percentages of A, T, G and C across the files read positions should be parallel. An imbalance in these equal proportions may suggest that there is an overrepresented sequence in your sample. 

**Per base GC content-** Provides a plot of the GC content (%) over each base position in the file. In theory, the overall base CG content would represent the GC content of the genome. In a random sample, the GC content would generally remain the same across the base pairs, so the line would run horizontally in this graph. A GC bias in a region of your base positions may suggest that your data contains an overrepresented sequence. An overall bias in GC content may suggest that there was a systematic issue during sequencing, or the sample library has a sequence bias. 

**Per sequence GC content-** Compares the GC content across the entire length of every sequence in the file to a modeled normal distribution of GC content for the genome. This model is created by using the observed data to create a reference genome. If you have a random library, one would expect the GC content to roughly follow a normal distribution, having the central peak corresponding to this modeled genome. A non-normal distribution shape could indicate contamination of the sample/ biased subset. A shift in distribution may suggest that there is a systematic bias occurring, that is independent of the base positions (*N.B: this error will not be flagged as orange/red by the program*).

**Per Base N Content-** When the sequencer cannot call a base with a high enough confidence, the base 'N' is called. This module represents the percentage of 'N' bases called for each base pair read position. 'N' calls often occur at the end of the read lengths; however, if the percentage of 'N' is greater than a few percent, you should question the validity of your base calls. 

**Sequence length distribution-** Produces a graph showing the distribution of read lengths in the analysed file. Some NGS sequencers generate uniform read length, whilst others do not. A warning will automatically be raised for this module if the reads in a file are of unequal length. 

**Duplicate sequences-** Highlights the level of sequence duplication in the sample file, by producing a plot showing the relative number of sequences with varying degrees of duplication. Note well, only the first 200, 000 reads in the file are analysed, and any sequence with more than 10 duplicates is automatically placed into the '10 duplicate' category. A low level of duplication may suggest a target sequence has high coverage, whereas a high level of duplication could indicate enrichment bias (eg. overamplification of PCR). 

**Overrepresented sequences-** This section lists any sequence that accounts for more than 0.1% of the total number of sequences in the library file. In a random library there is generally a diversity in sequences. If a sequence represents a higher than expected fraction of the whole (here greater than 0.1%), it may indicate the sequence is extremely biologically relevant, or that the library is either contaminated or not as random as it was thought to be. If a sequence is overrepresented, the program will try to match the sequence to a contaminants database. 

**Overrepresented Kmers-** Kmer analysis provides an indication of the levels of exactly repeated sequences within your sequence library. FasQC performs Kmer analysis by breaking up the first 20% of reads in your sequence library into 5-mers and then extrapolating the remaining portion of your library (fig.1). It calculates an observed to expected ratio for each Kmer, by first analysing the whole library to determine the expected level of each Kmer, and then comparing this to the observed Kmer counts. Any Kmer showing an overall 3 fold higher observed to expected ratio, or a 5 fold ratio at a specific base position, is reported by the module. This module provides a graph for the 6 Kmers with the highest number of hits, showing their pattern of enrichment across the length of your reads. This graph can be used to determine if you have general enrichment, or a pattern of bias across specific points within your read length (eg. due to adaptors). 

*To understand the concept of Kmers and how the Kmer counts are generated, refer to figure 1 below.* ![image label](https://github.com/BecWardell/NGS-pipelines/raw/bec_documentation/Docs/img/kmerCount.jpg)

> **Figure 1:** The mechanism by which FastQC analyses the Kmers in a sequence. The FastQC program analyses the first 20% of the library using this sliding 5-mer process. The remainder of the Kmer count is extrapolated based upon this analyses to produce a list of the highest Kmers in your library. 

**References:**

Babraham Bioinformatics. (2014) FastQC. Babraham Institute, Cambridgeshire, UK. Obtained from <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/> on the 06/04/2014. 
Wageningen University Bioinformatics Department. (2012) FastQC Manual. Wageningen University, Netherlands. Obtained from <http://www.bioinformatics.nl/courses/RNAseq/FastQC_Manual.pdf> on the 06/04/2014. 


Scythe
------

Function: 

Scythe uses a Naive Bayesian model to remove 3' adaptors (also called 'contaminants') from your reads, which you added when you labelled your samples prior to sequencing. Scythe works by looking at matches between sequences, and uses two probability models to determine if the match is due to contamination (adaptors) or by chance. The probability models are: i) the chance of their being these sequence matches given contamination or ii) the chance of these sequence matches being random. If the Bayesian probability model determines the read is contaminated due to high matching because it is likely an adaptor, that region of sequence will be removed. This model is possible because you have to set the rate of contamination in your sample (known as a prior), therefore this predicts the probability of your reads having matching sequences, and will consequently influence the Bayesian probability model. 


Use: 

To accurately remove 3' adaptors from each of your reads. Scythe is an accurate tool as it is able to account for the quality of your sequence, which generally decreases towards the 3' end where the adaptors are located. 


Requirements:

- You must set a prior (denoted as 'p'), which is your prediction of the adaptor contamination rate within your sample. In order to estimate your prior, and therefore set a limit for scythe, you can alter the script below.

```
    Example: given you adapter stars with ACCAGT
    adapt="ACCAGT"
    reads=100000
    fqfile=reads.fq
    echo "Contamination percent esitmate is: $(python -c print\ $(zcat $fqfile |head -n $reads| grep $adapt |wc -l)/${reads}.0*100)%"
```

- You must provide an adaptor sequence file so that the program knows what the contaminant sequences are. 


References:

UC Davis Genome Centre. (2014) Software. UC Davis Genome Center, Davis, California, USA. Obtained from <https://bioinformatics.ucdavis.edu/software/> on the 08/04/2014. 

Sickle
------

Function: 

Sickle acts to trim regions of reads that have deteriorating quality at the 3' and 5' ends (a result of NGS sequencing), by using set quality/length thresholds and a sliding windows technique. When each base is called by the sequencer, it is given a quality score. Sickle takes a specified window of bases (say window x, equal to 0.1 times read length) and finds the average quality score for this window. It then determines if the average quality score for window x is greater/less than a predetermined quality score threshold that you provide. If the average base calling quality is below the threshold, the sequence will be trimmed before the window (if at 3' end) or after the window (if at 5' end). In this way, Sickle determines when the base calling quality is sufficiently low to trim at the 3' end, and when the base calling quality has become high enough to trim at the 5' end of each read. If the sequence after trimming is below the specified read length threshold, the read will be discarded. If only one read within a pair is discarded, it is denoted as a 'single' read and the output placed into a seperate file. *Please refer to this figure to see a basic diagram of this window process.*


Use:

Most highthroughput NGS sequencing technologies produce reads with poor quality base calling at the 3' and sometimes 5' ends, which can be deleterious for downstream data analysis (assembly, mapping etc.). Thus, Sickle is used to remove poor quality 3' and 5' ends and reads that do not reach a specific length threshold.

*To show the difference in read quality before and after processing, you can run fastqc again*


Requirements:

Sickle requires you to set the following variables:

* -f: specify forward read file input.
* -r: specify reverse read file input.
* -q: the quality score limit (eg. if 20, the window must have an average quality score of 20 otherwise it will be discarded).
* -l: the minimum length threshold (eg. if 20, you would discard reads below 20 bps).
* -n: allows you to remove all sequences containing an n base. This feature is not recommended for RNAseq. 
* -t: the type of sequencer you used, either Sanger, Solexa or Illumina. This variable is extremely important because it is used to determine the quality scores of you bases. For Sanger, `33` is added to each quality score, whereas `64` is added to Solexa. Pre 2011, `33` is added to each quality score for Illumina, and post 2011 `64` is added.  
* -o: specify output file for forward reads. 
* -p: specify output file for reverse reads.
* -s: specify output file for non-paired reads. 


References:

UC Davis Genome Centre. (2014) Software. UC Davis Genome Center, Davis, California, USA. Obtained from <https://bioinformatics.ucdavis.edu/software/> on the 08/04/2014. 

Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files 
(Version 1.21) [Software].  Obtained from <https://github.com/najoshi/sickle> on the 08/04/2014. 


Subread
-------

Function: 

Subread uses a 'seed-and-vote' strategy to align reads back to a reference genome. This 'seed-and-vote' strategy essentially involves breaking up each read into several 'subreads' (16 nts long) and allowing each subread to vote on its optimum location in the genome. Thus, the region of the genome with the highest number of subreads theoretically corresponds to the region of the read (refer to fig.2 below). When reads are greater than 160 bp in length, overlapping subreads are used. Once the read location has been determined in the genome by having the highest number of subread 'votes', conventional algorithims fill in any in/del and mismatch information between the subreads. Thus, this tool is fast as the read is mapped onto the genome before any detailed filling in of missing regions in the read occurs. This tool is also sensitive, as it allows individual subreads to map to their optimum location meaning the subreads do not have to map close to eachother nor map exactly to the genome. It also has a high level of accuracy, requiring the final read location to correspond to several subreads. Subread however does have some limitations, for instance you cannot tell if alternative splicing has occurred, and therefore you should consider which genome mapping program is appropriate for your use. 


Use: 

To accurately and with great sensitivity align reads back to a reference genome. 


Requirements:

You must provide a reference genome that subread can use to map the reads onto. You need to use the subread built genome index, as subread relies upon the genome already been indexed. If you do not already have your genome indexed, please use the program 'subread-buildindex' so that subread can run. 

``` 
subread-buildindex -o <basename> -M <int> {FASTA FILE} [FASTA FILE 2]

#-o <basename>: name of index eg. TAIR10 
#int: how much memory the program uses, default 8GBs
#FASTA FILE: You must supply at least one file for the genome, additional files can be supplied using square brackets [file 2] etc. 
```

![image label](https://github.com/BecWardell/NGS-pipelines/raw/bec_documentation/Docs/img/seedAndVoteSubread.jpg)
> **Figure 2:** Artificial example of the subread ‘seed-and-vote’ mapping paradigm. (A) Here, six sub-reads are created from an artificial read each containing 5 continuous bases. The numbers in blue on the LHS denote where these sequences originated from in the read. The base sequence for each subread is encoded by a string of binary (seen in the hash table) as 2-bit binary is used to encode each base. (B) These subreads are then matched to perfectly complimentary regions in the reference genome (no mismatches allowed). Note that each subread may match to more than one location. Mapping identifies four candidate  locations, with 2, 5 1 and 2 votes respectively. The location with the greatest number of votes, in this case 5 votes, is chosen as the final mapping location for the read. (Modified from Liao Y et. al 2013). 

References:

Liao Y, Smyth GK, Shi W. (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. *Nucl. Acids Res*, 41(10):e108. 



FeatureCounts
--------------

Function:

FeatureCounts is a general purpose read-summarisation program for RNA and DNA sequencing, that identifies and quantifies the overall coverage for a specified genomic feature eg. genes, exons. FeatureCounts comparable in summarisation accuracy to other read counting programs, but has the major advantage of reduced computational cost. FeatureCounts has been found to be ten times faster on average than other current read counting methods and requires much less computational memory. This program is also the only one of its kind that can currently be run in parallel (eg. you can run all your files at once provided you have fewer files than cores on your computer). The high efficiency of featureCounts can be attributed to its ultrafast search algorithm for features, as well as the program being implemented in the C programming language. 


Use: 

To identify and then quantify reads for any genomic feature eg. determine expression level of genes, exons, transposons etc.  


Requirements:

 -You must provide the same genomic library used in subread for analysis, however, you need to provide the annotated genome version. 
 -The standard format for the annotated genome may be downloaded from library websites (eg. TAIR) in the standard gff3 format. You can however also use a saf file, which is a more basic condensed version of gff3. 


References:

Liao Y, Smyth GK, Shi W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923-930.

---

Miscellaneous
-------------

###Concatenate###

You can use this script when sequencing returns two forward reads and two reverse reads files for each sample. Essentially, the script merges the two forward reads files to create one forward reads file, and does likewise for the reverse reads. This produces single forward and reverse files for a sample that can then be modified in the rest of the pipeline prior to analysis. 

When used:

This script is used prior to the 01-runner.sh for fastqc. It is used when you have two forward reads files and two reverse reads files for a given sample. 

###BAM to TDF###

This script is used to convert a BAM file to a TDF file following featurecounts being run. It is useful as it reduces the size of the BAM file, making a smaller format file  for viewing gene expression levels in your genome browser.

When used: 

This scrip is used after featurecounts, which produces an output file in BAM format. 


###Results logs scraper###

Kevin's python script to extract key metrics out of the log files to make summary statistics. It is essential to check the metrics produced to ensure that the pipeline has correctly run. 
