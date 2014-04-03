#!/bin/sh
#  08-getScytheResults.sh
#
#
#  Created by Peter Crisp on 27/03/2014. Mauled to pieces by KDM.
#
set -e
set -x

logscript=$(readlink -f $0)/../tools/scythe_log.py

find logs/reads_noadapt.*/  -mindepth 1 -maxdepth 1 -type f -name Sample* | xargs python $logscript >results_summaries/reads_noadapt.csv
