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
06-runner.sh <strandedness>"

######### Setup ################
strand=$1
alignFolder=$2
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "2" ]
then
echo $usage
exit -1
else
echo "featureCounts strandedness setting = $1\n alignment folder = $2"
fi
########## Run #################


#user defined variables that could be changed:
workingdir=./
script=$scriptdir/06-featureCounts.sh
###

function findSamples () {
find ${alignFolder}/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=featureCounts
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel bash $script {} $strand $alignFolder \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#usage:
#bash ~/path_to/06-runner.sh <strandedness # >
#for featureCounts three values are possible: o =unstranded, 1= stranded, 2= reversely stranded.

