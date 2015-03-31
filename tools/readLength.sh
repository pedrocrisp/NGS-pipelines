#!/bin/sh

#  readLength.sh
#  get read lengths of all reads ina fastq and make a text histogram
#
#  Created by Peter Crisp on 28/03/2015.
#


#Using perl:
#
#cat input.fq | perl -ne '$s=<>;<>;<>;chomp($s);print length($s)."\n";' > input.readslength.txt
#
#Using awk:
#
#cat input.fq | awk '{if(NR%4==2) print length($1)}' > input.readslength.txt
#
#if zipped file, using:
#
#zcat input.fq.gz | ...
#
#get length statistics:
#
#sort input.readslength.txt | uniq -c
#
#or using fancy code from Jim Kent's utility:
#
#textHistogram
#
#So, one line code for all input fastq files would be:

##########

## to run seperately and save the output:

# zcat *.*q.gz | awk '{if(NR%4==2) print length($1)}' > readslength.txt
# textHistogram -maxBinCount=59 readslength.txt > readslength.summary.txt


find *.fq.gz -not -name \*raw\* -printf "zcat %p | awk '{if(NR%%4==2) print length(\$1)}' | textHistogram -maxBinCount=59 stdin \n" | sh