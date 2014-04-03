#!/bin/sh

#  08-getfeatureCountsResults.sh
#  
#
#  Created by Peter Crisp on 27/03/2014.
#

set -e
set -x

logscript=$(dirname $(readlink -f $0))/../tools/featureCounts_log.py
find logs/featureCounts.*/  -mindepth 1 -maxdepth 1 -type f -name Sample* | xargs python $logscript >results_summaries/featureCounts.csv
