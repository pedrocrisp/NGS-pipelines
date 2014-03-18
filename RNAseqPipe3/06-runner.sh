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

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/06-featureCounts.sh
###

function findSamples () {
find align/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=featureCounts
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel bash $script {} \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/06-runner.sh
