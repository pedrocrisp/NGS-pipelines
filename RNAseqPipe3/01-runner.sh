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
script=$scriptdir/01-fastqc.sh
outdir=reads_fastqc
###

function findSamples () {
find reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

mkdir logs
logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel bash $script {} \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run, got to directory containing reads directory and call:
#bash ~/path_to/01-runner.sh
