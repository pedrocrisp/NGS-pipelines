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
05-runner.sh <number of threads> <bam folder> <gff reference> <genelist to subset gff>"

######### Setup ################
threads=$1
bam_folder=$2
gff_ref=$3
genelist=$4

if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel bowtie jobs on $reads folder, bowtie2 can use $b_threads threads"
fi
########## Run #################

script=$scriptdir/05-bowtie2.sh
###

function findSamples () {
find $bam_folder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

outdir=${bam_folder}_parestahp
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel -j $threads bash $script {} $bam_folder $gff_ref $genelist $outdir \>logs/${outdir}_bowtie2.${timestamp}/{}.log 2\>\&1

#To run:
#bash ~/path_to/08-runner.sh
