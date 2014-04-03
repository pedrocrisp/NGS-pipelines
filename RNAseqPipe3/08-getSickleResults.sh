#!/bin/sh

#  08-getSickleResults.sh
#
#
#  Created by Peter Crisp on 27/03/2014.
#

set -e
set -x

logscript=$(dirname $(readlink -f $0))/../tools/sickle_log.py
find logs/reads_noadapt_trimmed.*/  -mindepth 1 -maxdepth 1 -type f -name Sample* | xargs python $logscript >results_summaries/sickle.csv
