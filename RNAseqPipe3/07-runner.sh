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
06-runner.sh <strandedness> <alignment folder>"

######### Setup ################
alignFolder=$1
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "1" ]
then
echo $usage
exit -1
else
echo "alignment folder = $1"
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

findSamples | parallel bash $script {} $alignFolder \>logs/${outdir}_subread.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/07-runner.sh
