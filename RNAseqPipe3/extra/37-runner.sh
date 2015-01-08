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
37-runner.sh <exon_beds folder> <threads> <min mRNA Length (nt)> <trim Length (nt)> <min mRNA coverage>
suggested min mRNA length 600nt
suggested trim length 250 (from each end)
suggested min average per bp depth of coverage 100?"

######### Setup ################
alignFolder=$1
threads=$2
minLength=$3
trimLength=$4
minCoverage=$5

if [ "$#" -lt "5" ]
then
echo $usage
exit -1
else
echo "making exon bd files \n alignment folder = $1\n iniating $2 parallel 27-bam_to_exon_stranded jobs\n minimum transcript length = $3\n timming $4 bases from 3' and 5' ends\n minimum average per base coverage $5"
fi
########## Run #################


#user defined variables that could be changed:
workingdir=./
script=$scriptdir/37-mRNA-density.R
###

function findSamples () {
find ${alignFolder}/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir="${alignFolder}_plots"
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads R -f $script --args {} $alignFolder $outdir $minLength $trimLength $minCoverage \>${logdir}/{}.log 2\>\&1

#To run:
#Must be run after 27 which generates the .bed.gz files
#bash ~/path_to/37-runner.sh
