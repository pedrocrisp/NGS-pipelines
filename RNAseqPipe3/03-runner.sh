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
03-runner.sh <encoding> <number of threads> <SE or PE>"

######### Setup ################
encoding=$1
threads=$2
read_ends=$3
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel quality trimming jobs using sickle, $2 encoding, or $3 reads"
fi
########## Run #################

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/03-sickle.sh
outdir=reads_noadapt_trimmed
###

function findSamples () {
find reads_noadapt/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

mkdir $outdir
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $encoding $read_ends \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/03-runner.sh
