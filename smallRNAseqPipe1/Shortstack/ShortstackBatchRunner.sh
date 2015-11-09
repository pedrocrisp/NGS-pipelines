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
03-runner.sh <number of threads to parallel> <fastq folder>"

######### Setup ################
threads=$1
fastqFolder=$2

if [ "$#" -lt "2" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel Shortstack jobs"
fi
########## Run #################

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/ShortstackBatch.sh
outdir=shortstacked
###

function findSamples () {
find $fastqFolder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

mkdir $outdir
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $fastqFolder $outdir \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/03-runner.sh
