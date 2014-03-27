#!/bin/sh

#  08-getSubreadResults.sh
#  
#
#  Created by Peter Crisp on 27/03/2014.
#

set -e
set -x

function findSamples () {
find logs/align_subread.*/  -mindepth 1 -maxdepth 1 -type f -name Sample* | tr ' ' '\n'
}

samples=$(findSamples)

for sample in $samples
do
echo $(basename $sample)
grep -e Processed -e Mapped -e paired -e Indels $sample
done
