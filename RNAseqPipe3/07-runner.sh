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
07-runner.sh <alignment folder> <number of threads> <strandedness of library>"

######### Setup ################
alignFolder=$1
threads=$2
strand=$3
if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "alignment folder = $1, initiating $2 parallel make stranded bigWig jobs, library is $3"
fi
########## Run #################

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/07-bam_to_tdf_stranded.sh
###

function findSamples () {
find $alignFolder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=tdf_for_igv
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}_subread.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $alignFolder $strand \>logs/${outdir}_subread.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/07-runner.sh
