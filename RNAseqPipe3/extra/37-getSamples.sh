#!/bin/sh

#  37-getSamples.sh
#  
#
#  Created by Peter Crisp on 7/10/2014.
#

set -e
set -x

sample=$1
sample_dir=exon_beds/$sample
cd $sample_dir

R -f 