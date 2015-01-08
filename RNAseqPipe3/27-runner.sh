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
27-runner.sh <alignment folder> <threads>"

######### Setup ################
alignFolder=$1
threads=$2

if [ "$#" -lt "2" ]
then
echo $usage
exit -1
else
echo "making exon bd files \n alignment folder = $2\n iniating $3 parallel 27-bam_to_exon_stranded jobs"
fi
########## Run #################

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/27-bam_to_exon_stranded.sh
###

function findSamples () {
find ${alignFolder}/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

resultsFolder="${alignFolder}_exon_beds"
mkdir $resultsFolder
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${resultsFolder}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $alignFolder $resultsFolder \>${logdir}/{}.log 2\>\&1

#To run:
#Must be run after 07 so that the split strand bams are created
#bash ~/path_to/27-runner.sh <alignment folder> <threads>
