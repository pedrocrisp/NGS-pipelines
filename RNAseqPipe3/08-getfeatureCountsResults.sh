#!/bin/sh

#  08-getfeatureCountsResults.sh
#  
#
#  Created by Peter Crisp on 27/03/2014.
#

set -e
set -x

function findSamples () {
find logs/featureCounts.*/  -mindepth 1 -maxdepth 1 -type f -name Sample* | tr ' ' '\n'
}

samples=$(findSamples)

for sample in $samples
do
echo $(basename $sample)
grep 'Number of successfully assigned fragments is' $sample
done