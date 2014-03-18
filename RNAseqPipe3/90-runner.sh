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
script=$scriptdir/90-head2.sh
#outdir=~/workspace/test_data/reads
###

function findSamples () {
find reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

#mkdir ${outdir}

findSamples | parallel bash $script {}

#To run, got to directory containing reads directory and call:
#note original full sized files will be deleted
#bash ~/path_to/90-runner.sh
