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
script=$scriptdir/05-subread.sh
###

function findSamples () {
find reads_scythe_seqtk/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=align
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}_subread.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel bash $script {} \>logs/${outdir}_subread.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/05-runner.sh
