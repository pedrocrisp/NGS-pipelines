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
GeneCoverageR-runner.R <distance beds folder> <threads>"

######### Setup ################
alignFolder=$1
threads=$2
library_layout=$3
#minLength=$3
#trimLength=$4
#minCoverage=$5

if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "making Gene Coverage plots \n alignment folder = $1\n iniating $2 parallel jobs"
fi
########## Run #################


#user defined variables that could be changed:
workingdir=./
script=$scriptdir/GeneCoverageR.R
###

function findSamples () {
    find $alignFolder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir="${alignFolder}_plotData"
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="${outdir}/.logs_${timestamp}"
mkdir $logdir

cat $script > "${logdir}/script.log"
cat $0 > "${logdir}/runner.log"
cat $script

findSamples | parallel -j $threads R -f $script --args {} $alignFolder $outdir $library_layout \>${logdir}/{}.log 2\>\&1

#To run:
#bash ~/gitrepos/NGS-pipelines/GeneCoverage/GeneCoverageR-runner.sh tdf_for_igv_coverage_beds 1
