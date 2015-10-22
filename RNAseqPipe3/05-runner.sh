#!/bin/bash

#  05-runner.sh
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
Tflag=
aflag=
iflag=
Pflag=
while getopts j:T:a:i:P: name
do
    case $name in
        j)    jflag=1
        jval="$OPTARG";;
        T)    Pflag=1
        Tval="$OPTARG";;
        a)    aflag=1
        aval="$OPTARG";;
        i)    iflag=1
        ival="$OPTARG";;
        P)    Pflag=1
        Pval="$OPTARG";;
        ?)   printf "Usage: %s: [-j value] [-T value] [-a <aligner>] [-i <reference>] [-P value]\n" $0
            exit 1;;
        *)   printf "Usage: %s: [-j value] [-T value] [-a <aligner>] [-i <reference>] [-P value]\n" $0
            exit 1;;
    esac
done
if [ ! -z "$jflag" ]; then
    printf 'Parallel will use -j "%s" threads\n' "$jval"
fi
if [ ! -z "$Tflag" ]; then
    printf 'subread will use -T "%s" threads\n' "$Pval"
fi
if [ ! -z "$aflag" ]; then
    printf '"%s" will be used for alignment\n' "$aval"
fi
if [ ! -z "$iflag" ]; then
    printf '"%s" will be used for reference index\n' "$ival"
fi
if [ ! -z "$Pflag" ]; then
printf 'subread will assume -P "%s" encoding\n' "$Pval"
fi
shift $(($OPTIND - 1))

if [ -z "$jflag" ]
then
printf "Threads for parallel not specified\n Usage: %s: [-j value] [-T value] [-a <aligner>] [-i <reference>] [-P value]\n" $0
exit
fi

if [ -z "$Pflag" ]
then
printf "Threads for subread not specified\n Usage: %s: [-j value] [-T value] [-a <aligner>] [-i <reference>] [-P value]\n" $0
exit
fi

if [  "$aval" == "subread-align" ]
then
printf "aligner subread-align\n"
elif [ "$aval" == "subjunc" ]
then
printf "aigner subjunc\n"
else
printf "Aligner not specified\n Usage: %s: [-j value] [-P value] [-a <aligner>]\n" $0
exit
fi

if [ -z "$iflag" ]
then
printf "Refernec index for subread not specified\n Usage: %s: [-j value] [-P value] [-a <aligner>] [-i <reference>]\n" $0
exit
fi

if [ -z "$Pflag" ]
then
printf "Encoding of fastQ not specified\n Usage: %s: [-j value] [-T value] [-a <aligner>] [-i <reference>] [-P value]\n" $0
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
script=$scriptdir/05-subread.sh
###

function findSamples () {
find reads_noadapt_trimmed/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir="$aval"
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $jval bash $script {} $Tval $aval $ival $Pval \>$logdir/{}.log 2\>\&1

#To run:
#bash ~/path_to/05-runner.sh
