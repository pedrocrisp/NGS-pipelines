RNAseq pipeline user guide
===========================

Fast-QC
-------

- Quality control
- Generates a http file
- Want to ensure that your reads are mostly in the green section.
- Why do? Because your read quality decreases at the end of each read.
- Here is an example of a [good quality](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc/fastqc_report.html) read, note the 'per base sequence quality' has very low variation and across the reads it has a high quality score. You should look at each of the features for a 'tick' as these vary depending on your read quality. 


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
