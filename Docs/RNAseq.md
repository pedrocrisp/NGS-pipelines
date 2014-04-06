RNAseq pipeline user guide
===========================

Fast-QC
-------

Function: 

Quality control of raw/altered data from high throughput sequencing. 
FastQC is a java based program designed to perform a quality control analysis on your reads to identify any issues in the sequencing tool (eg. read quality decreases at the end of every read) / or in your starting library material. 

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

Per base sequence quality: Shows the variation in sequence quality at each base position for all reads in the file. The y axis represents the quality scores and the x axis is the base position in the reads. The information is presented as a boxplot for each read position, with the red line representing the median and the blue line representing the mean values. As quality of the bases increases with quality number, the graph has been divided into 
Per sequence quality scores:
Per base sequence content:
Per base GC content:
Per sequence GC content:


read, note the 'per base sequence quality' has very low variation and across the reads it has a high quality score. 

References:

Babraham Bioinformatics. (2014) FastQC. Babraham Institute, Cambridgeshire, UK. Obtained from <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/> on the 06/04/2014. 
Wageningen University Bioinformatics Department. (2012) FastQC Manual. Wageningen University, Netherlands. Obtained from <http://www.bioinformatics.nl/courses/RNAseq/FastQC_Manual.pdf> on the 06/04/2014. 




Scythe
------

- What? Removes adaptors
- How? Does a simple match, and if it finds mis-matches it uses a basian calculation to calculate the probability of a match to an adaptor and trims if it is a match. 
- Have to set a prior (aka: prediction of the contamination rate). Probability of it matching influences how many it says will match. Aim to set to what you think. Calculate vague contamination ratio. Adaptor is contamination.  
-

Sickle
------

*To show the difference in read quality, fastqc can be run again*

Subread
-------

Feature counts
--------------

---

Miscellaneous
-------------

###Concatenate###

###BAM to TDF###

###Results logs scraper###
