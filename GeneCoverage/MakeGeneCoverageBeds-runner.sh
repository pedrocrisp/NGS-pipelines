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
MakeGeneCoverageBeds-runner.sh <alignment folder> <number of threads> <reference>
Reference might be ~/ws/refseqs/TAIR10/TAIR10_crisp_gene_primary.bed
"

######### Setup ################
alignFolder=$1
threads=$2
reference=$3
if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "alignment folder = $1, initiating $2 parallel MakeGeneCoverageBeds jobs, reference is $3"
fi
########## Run #################

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/MakeGeneCoverageBeds.sh
###

function findSamples () {
find $alignFolder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir="${alignFolder}_coverage_beds"
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./${outdir}/.logs_${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $alignFolder $reference $outdir \>${logdir}/{}.log 2\>\&1

#To run:
#bash ~/path_to/MakeGeneCoverageBeds-runner.sh
