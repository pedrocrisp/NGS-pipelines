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
#

usage="USAGE:
05-runner.sh <number of threads> <reads folder> <bowtie threads per job> <reference> <number of mapping locations to report>"

######### Setup ################
threads=$1
reads=$2
b_threads=$3
reference=$4
multimapping=$5
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "5" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel bowtie jobs on $reads folder, bowtie2 can use $b_threads threads"
fi
########## Run #################

script=$scriptdir/05-bowtie2.sh
###

function findSamples () {
find $reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=align_bowtie2
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}_bowtie2.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $reads $b_threads $reference $multimapping \>logs/${outdir}_bowtie2.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/05-runner.sh
