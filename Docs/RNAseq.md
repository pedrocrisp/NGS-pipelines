RNAseq pipeline user guide
===========================

Fast-QC
-------

Function: 

FastQC is a java based program designed to perform a quality control analysis on your reads to identify any issues that occurred during high throughput sequencing (eg. read quality decreases at the end of every read) / or in your starting library material. 

Use: 

Fast QC can be used to perform a simple quality control analysis of your raw data, and also be used after modifications (eg. removal of adapter sequences) to see if any improvement in quality has occurred. 

Requirements:

Fast QC can read files ending in: FastQ (all quality encoding variants), Casava FastQ files^, Colorspace FastQ, GZip compressed FastQ, SAM, BAM or SAM/BAM Mapped only^. 
^Indicates that you will need to tell FastQC you are using SAM/BAM or Casava files, as FastQC will not select these automatically. 

Evaluation of results:

- Once you have run Fast-QC, the results report will be available as a http file. When you click on it, the report should open in your browser.
- You will see that the analysis is segmented into a series of modules, with a summary of the module outcomes provided on the left hand side. Ratings include: normal (green tick), slightly abnormal (orange exclamation mark) and very abnormal (red cross). 
- Here is an example of a [good quality](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc/fastqc_report.html) and [poor quality] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc/fastqc_report.html) read. A brief description of the purpose of each module is provided below. 

Basic statistics module: Provides a summary of the data, including filename/type, total sequences (actual or estimate??? ASK KEV), sequence length (provides the range of read lengths, if all reads are the same length, only one number will appear) and overall GC content expressed as a percentage. N.b: basic statistics will always appear to have a green tick (never flags a warning/error), so you should always further analyse the data. 

Per base sequence quality: Shows the variation in sequence quality at each base position for all reads in the file. The y axis represents the quality scores and the x axis is the base position in the reads. The information is presented as a boxplot for each read position, with the red line representing the median and the blue line representing the mean values. As base calling quality increases with quality number, the graph has been divided into green, orange and red sections, respectively representing good, satisfactory and poor quality base calls. In a good quality run, one would want their boxplots to be in the green section and have little variation. 

Per sequence quality scores: enables you to see if a subset of your reads have overall poor quality base recall. Often a small percentage of your reads will be universally poor quality, generally because they have been poorly imaged (eg. may be on edge of the field of view). If there is a significant percentage of reads that have poor quality in a region of base pairs, this may suggest a systematic problem. 

Per base sequence content: provides a plot of the percantage of each base called (A, T, G, C) at each base pair position in the file. In theory, in a random library, the proportion of the four bases for a given base pair position should be approximately equal, thus the lines for A, T, G, C percentages should be in parallel. An imbalance in these equal proportions may suggest that there is an overrepresented sequence in your sample. 

Per base GC content: provides a plot of the GC content (%) over each base position in the file. In a random sample, the GC content would generally remain the same across the base pairs (so run horizontally). In theory, the overall CG content would represent the GC content of the genome. A GC bias in a region of your base positions may suggest your data contains an overrepresented sequence. An overall bias in GC content may suggest that there was a systematic issue during sequencing, or the sample library has a sequence bias. 

Per sequence GC content: compares the GC content across the entire length of every sequence in the file to a modeled normal distribution of GC content for the genome (model is created by using the observed data to create a reference genome). If you have a random library, one would expect the GC content to roughly follow a normal distribution, having the central peak corresponding to this modeled genome. A non-normal distribution shape could indicate contamination of the sample/ biased subset, and a shifted distribution could suggest that there's a systematic bias independent of the base positions (the latte will not be flagged as orange/red by the program).

Per Base N Content: When the sequencer cannot call a base with a high enough confidence, the base 'N' is called. This module represents the percentage of 'N' bases called for each base pair read position. 'N' calls often occur at the end of the read lengths, however, if the percentage of N is greater than a few percent, you should question the validity of your base calls. 

Sequence length distribution: produces a graph showing the distribution in read lengths in the analysed file. Some NGS sequencers generate uniform read length, whilst others do not. A warning will automatically be raised for this module if the reads in a file are of unequal length. 

Duplicate sequences: Highlights the level of sequence duplication in the sample file, by producing a plot showing the relative number of sequences that have varying degrees of duplication. Note well, any sequence with more than 10 duplicates is automatically placed into the '10 duplicate' category and only the first 200, 000 reads in the file are analysed. A low level of duplication may suggest a sequence has high coverage, whereas a high level of duplication could indicate a sequence has been enriched eg. overamplification of PCR. 

Overrepresented sequences: This section lists any sequence that accounts for more than 0.1% of the total number of sequences. In a random library there is generally a diversity in sequences. If a sequence represents a higher than expected fraction of the whole (here greater than 0.1%), it may suggest that the sequence is extremely biologically relevant, the library is contaminated, or not as random as thought. If a sequence is overrepresented, the program will try to match the sequence to a contaminants database. 

Overrepresented Kmers: Kmer analysis provides an indication of the levels of exactly repeated sequences within your sequence library. FasQC performs Kmer analysis by breaking up the first 20% of reads in your sequence library into 5-mers and then extrapolating the remaining portion of your library. It calculates an observed/expected ratio for each Kmer by determinng the expected level for each Kmer based upon the base content of the whole library and comparing this to the observed Kmer counts. Any Kmer showing an overall 3 fold observed/expected ratio or a 5 fold ratio at a specific base position is reported by the module (these are the threshold conditions). It also draws a graph for the 6 Kmers with the highest number of hits, showing their pattern of enrichment across the length of your reads. This graph can be used to determine if you have general enrichment or bias within your read length (eg. due to adaptors)

*To understand the concept of Kmers and how the Kmer counts are generated, refer to this figure.* *INSERT FIG*

References:

Babraham Bioinformatics. (2014) FastQC. Babraham Institute, Cambridgeshire, UK. Obtained from <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/> on the 06/04/2014. 
Wageningen University Bioinformatics Department. (2012) FastQC Manual. Wageningen University, Netherlands. Obtained from <http://www.bioinformatics.nl/courses/RNAseq/FastQC_Manual.pdf> on the 06/04/2014. 


Scythe
------

Function: 

Scythe uses a Naive Bayesian model to remove 3' adaptors (also called 'contaminants') from your reads, which you added when you labelled your samples prior to sequencing. Scythe works by looking at matches between sequences, and uses two probability models to determine if the match is due to contamination (adaptors) or by chance. The probability models are: i) the chance of their being these sequence matches given contamination or ii) the chance of these sequence matches being random. If the Bayesian probability model determines the read is contaminated due to high matching because it is likely an adaptor, that region of sequence will be removed. This model is possible because you have to set the rate of contamination in your sample (known as a prior), therefore this predicts the probability of your reads having matching sequences, and will consequently influence the Bayesian probability model. 


Use: 

To accurately remove 3' adaptors from each of your reads. Scythe is an accurate tool as it is able to account for the quality of your sequence, which generally decreases towards the 3' end where the adaptors are located. 

Requirements:

You must set a prior (denoted as 'p'), which is your prediction of the adaptor contamination rate within your sample. In order to estimate your prior, you TEXT HELP TEXT HELP. 

    # say your adapter stars with ACCAGT
    adapt="ACCAGT"
    reads=100000
    fqfile=reads.fq
    echo "Contamination percent esitmate is: $(python -c print\ $(zcat $fqfile |head -n $reads| grep $adapt |wc -l)/${reads}.0*100)%"


References:

https://bioinformatics.ucdavis.edu/software/

Sickle
------

Function: 

Sickle acts to trim regions of reads that have deteriorating quality at the 3' and 5' ends (a result of NGS sequencing), by using set quality/length thresholds and a sliding windows technique. Sickle determines when the base calling quality is sufficiently low to trim at the 3' end, and when the base calling quality has become high enough to trim at the 5' end of each read. Sickle also trims reads based upon read length, FINISH THIS AND DRAW A DIAGRAM!


Don't comprehened:
"It takes the quality values and slides a window across them whose length is 0.1 times the length of the read. If this length is less than 1, then the window is set to be equal to the length of the read. Otherwise, the window slides along the quality values until the average quality in the window rises above the threshold, at which point the algorithm determines where within the window the rise occurs and cuts the read and quality there for the 5′-end cut. Then when the average quality in the window drops below the threshold, the algorithm determines where in the window the drop occurs and cuts both the read and quality strings there for the 3′-end cut. However, if the length of the remaining sequence is less than the minimum length threshold, then the read is discarded entirely."

Use:

Most highthroughput NGS sequencing technologies produce reads with poor quality base calling at the 3' and sometimes 5' ends, which can be deleterious for downstream data analysis (assembly, mapping etc.). Thus, Sickle is used to remove poor quality 3' and 5' ends and reads that do not reach a specific length threshold.

Requirements:

What do I require? 
q= quality score limit
l= length score limit
n= you can remove all sequences with n - wouldn't recommend for RNAseq. 
t= sanger quality score and coding (33+). Solexa (64+), Illumina (pre 2011 33+, after 64+)

References:
https://bioinformatics.ucdavis.edu/software/

*To show the difference in read quality, fastqc can be run again*

Subread
-------

Function: Subread uses a 'seed-and-vote' strategy to align reads back to a reference genome. This 'seed-and-vote' strategy essentially involves breaking up each read into several 'subreads' (hence the name) and allowing each subread to vote on its optimum location in the genome. Thus, the region of the genome with the highest number of subreads theoretically corresponds to the region of the read. When reads are greater than 160 bp in length, overlapping subreads are used. Between the subreads that won the vote and have therefore determined the read location on the genome, conventional algorithims fill in any in/del and mismatch information. This tool is fast as mapping of the read onto the genome is performed prior to filling in detail of missing reads. Being flexible to let individual subreads map to their optimum location, and then determine the greatest voting block also means it is sensitive, as it does not require the subreads to map close to eachother or subreads to map exactly to the genome. This tool is also accurate, as the final location for the read must correspond to several subreads. 

Use: to accurately and with great sensitivity align reads back to a reference genome. 

Requirements:

References:
The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote

Feature counts
--------------

Function:

Use: A read-summarisation program to identify and then quantify reads for any genomic feature eg. exons.  

Requirements:

References:
http://bioinformatics.oxfordjournals.org/content/30/7/923.long

---

Miscellaneous
-------------

###Concatenate###

###BAM to TDF###

###Results logs scraper###
