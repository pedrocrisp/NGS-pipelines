#!/bin/sh

#  bam_lengths.sh
#  
#
#  Created by Peter Crisp on 25/05/2014.
#

if [ $# -lt 1 ]
then
echo "Usage: $0 [directory with bams] ..."
exit 1
fi

dir=$1

function findBams () {
find $dir -type f -name "*.bam" -exec echo {} \;| tr ' ' '\n'
}

echo "$0 counts the lines in the bam"
l=0
n=0
r=0

for bam in $(findBams)
do
l=`samtools view $bam | wc -l`
echo "$bam: $l"
n=$[ $n + 1 ]
r=$[ $r + $l ]
done

echo "$n files in total, with $r lines in total"
