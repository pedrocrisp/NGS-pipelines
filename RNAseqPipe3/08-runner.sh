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
outdir=results_summaries
###

mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $0 > "$logdir/$(basename $0)"

#scythe
getScytheResults=$scriptdir/08-getScytheResults.sh
cat $getScytheResults > "${logdir}/$(basename $getScytheResults)"
bash $getScytheResults > ${outdir}/scytheResults.txt

#sickle
getSickleResults=$scriptdir/08-getSickleResults.sh
cat $getSickleResults > "${logdir}/$(basename $getSickleResults)"
bash $getSickleResults > ${outdir}/sickleResults.txt

#subread
getSubreadResults=$scriptdir/08-getSubreadResults.sh
cat $getSubreadResults > "${logdir}/$(basename $getSubreadResults)"
bash $getSubreadResults > ${outdir}/subreadResults.txt

#featureCounts
getfeatureCountsResults=$scriptdir/08-getfeatureCountsResults.sh
cat $getfeatureCountsResults > "${logdir}/$(basename $getfeatureCountsResults)"
bash $getfeatureCountsResults > ${outdir}/featureCountsResults.txt


#findSamples | parallel bash $script {} \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run, got to directory containing reads directory and call:
#bash ~/path_to/08-runner.sh
