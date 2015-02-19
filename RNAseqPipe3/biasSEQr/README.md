# biasSEQr
##The Degradome Profiler pipeline
The Degradome Profiler pipeline interrogates RNA sequencing data (mRNAseq or PARE, although it could also be applied to any sequence data, for instance (epi)genomic sequencing or CHIPseq) for signatures of decay intermediates. The workflow of the pipeline is outlined below. In steps 1-4 reads are aligned to the genome and summarized into bed files (these steps could be customised to employ tools of preference). Here, reads are aligned to the Arabidopsis genome using *subjunc* from the subread package (*Liao, 2013*). It is critical to use a splice junction aware mapper to align the entire read or else artefacts - namely low coverage islands - will be produced at exon boundaries (for instance *subread-align* will trim reads in this way). In Step 2, the alignment files are split into reads mapping to the Watson and Crick strands in order to preserve the strand specificity (at the time of analysis bedtools *coverageBed* could not properly interpret paired-end strand specific data).   In Step 3, reference bed files for each strand with the co-ordinates of each exon (or other gene feature) are prepared from a reference GFF using *gff2bed* (bedops, requires python). For the present analysis only exons were considered. In Step 4, bed files summarising reads to coverage per base across the transcriptome are generated using *coverageBed* (bedtools). One file for each strand is created then both are concatenated together and compressed. 

Example work flow using the scripts in the RNAseqPipe3 directory is:

1. 01-runner.sh FASTQC
2. 02-runner.sh adapter removal scythe
3. 03-runner.sh quality trimming sickle
4. 04-runner.sh FASTQC
2. 05-runner.sh alignment sunjunc
2. 05b-runner.sh sam to bam (indexed and sorted)
2. 27-runner.sh strand specific bams
3. biasSEQr-runner.sh biasSEQr

![Degradome Profiler Pipeline](https://github.com/pedrocrisp/NGS-pipelines/blob/master/RNAseqPipe3/biasSEQr/DegradomeProfilerPipe.png)

## biasSEQr
The bed files produced by the front-end wrappers of the Degradome Profiler Pipeline are then analysed by the biasSEQr software in R to interrogate the distribution of reads along the length of all genes in the transcriptome. Reads are filtered by length, the ends trimmed to account for the reduced ability to capture the ends of a transcript during library preparation, and then further filtered by average per base coverage. Each transcript is then divided into 2, 10 and 100 fixed width windows to cover the length of a transcript. Base pair coverage per bin is summarised using *stats.bin()* and expressed as a percentage of the total base-pair counts per transcript yielding normalised coverage profiles. biasSEQr outputs heatmaps, coverage plots and summary csv files. Suggested setting for analysis, transcripts 401 nt or longer were trimmed 150 nt from each end and only those with an average base-pair coverage of 10 were retained for analysis. 

