#!bin/bash

#  cat_fastqs-runner.sh
#
#
#  Created by Peter Crisp on 26/3/2015.
#

set -e
set -x


#For use on multi fastq files ie L001 and L002 where sequences files are provided as multiple files NOT R2 and R2 which would be paired end data

###
#code to make it work on osx and linux

if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi

###

usage="USAGE:
cat_fastqs-runner.sh <number of threads> <reads folder> <PE or SE>"

######### Setup ################
threads=$1
reads=$2
PEorSE=$3
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel cat jobs on reads in $reads folder on $PEorSE data"
fi
########## Run #################

#cat_fastqs.sh should be in same folder as the runner:
script=$scriptdir/cat_fastqs.sh
###

function findSamples () {
find $reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

findSamples | parallel -j $threads bash $script {} $reads $PEorSE

#To run:
#bash ~/path_to/cat_fastqs-runner.sh
