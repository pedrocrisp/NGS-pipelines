#!/bin/bash

#  05b-runner.sh
#  subread runner script
#
#  Created by Peter Crisp on 14/10/2014.
#

#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced.
set -x

#getops
unset name
jflag=
Fflag=
while getopts j:F: name
do
    case $name in
        j)    jflag=1
        jval="$OPTARG";;
        F)    Fflag=1
        Fval="$OPTARG";;
        ?)   printf "Usage: %s: [-j value] [-F folder]\n" $0
        exit 1;;
        *)   printf "Usage: %s: [-j value] [-F folder]\n" $0
        exit 1;;
    esac
done
if [ ! -z "$jflag" ]; then
    printf 'Parallel will use -j "%s" threads\n' "$jval"
fi
if [ ! -z "$Fflag" ]; then
printf 'Folder with .sam files is "%s"\n' "$Fval"
fi
shift $(($OPTIND - 1))

if [ -z "$jflag" ]
then
    printf "Threads for parallel not specified\n Usage: %s: [-j value] [-F folder]\n" $0
exit
fi

if [ -z "$Fflag" ]
then
printf "Folder with .sam files not specified\n Usage: %s: [-j value] [-F folder]\n" $0
exit
fi



###find script directory (for purpose of locating target scripts and reference sequences
#code to make it work on osx and linux
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#set script dir and working dir
workingdir=./
script=$scriptdir/05b-samtools.sh
###

function findSamples () {
find $Fval/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=makeIndexedBams
#mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $jval bash $script {} $Fval \>$logdir/{}.log 2\>\&1

#To run:
#bash ~/path_to/05b-runner.sh
