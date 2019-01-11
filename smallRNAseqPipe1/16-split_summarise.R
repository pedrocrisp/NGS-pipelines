#!/usr/bin/Rscript
##########

# Peter Crisp
# 2018-08-06
# multimapping summary stats and filter
# this step creates an NH:i tag in R, summarises NH:i distribution and then removes reads with NH:i equal to multi_filter

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
multi_filter <- args[3]
multi_filter
results_folder <- args[4]
results_folder

###########################
# debug
sample="B73_S"
data_folder="tmp_out_dir"
multi_filter=11
results_folder="split_align"

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

# There are probably better ways to parse amd filter a sam but oh well
# I am be horribly mangling the sam format...
# this probably looses the headders - consider copying from the original files

##### 21 nt
sam_in <- read_tsv(paste0(data_folder, "/", sample, "/",  sample, "_q10_21.sam"), col_names = F)
# calculate number of alignments per read
out_sam <- sam_in %>% select(1:11) %>% group_by(X1) %>% mutate(MM = paste0("NH:i:", n()), n = n())

# make summary file
out_sam_summary <- out_sam %>% ungroup() %>% distinct(X1, .keep_all = T) %>% group_by(n) %>% summarise(multi = n())
out_sam_summary
write_tsv(out_sam_summary, paste0(results_folder, "/", sample, "_q10_21_summary.tsv"))

# filter out reads multimapping over the specified limit (multi_filter)
out_sam <- out_sam %>% select(-n) %>%
write_tsv(out_sam, paste0(results_folder, "/", sample, "/", sample, "_q10_21_filtered.sam"), , col_names=F)

##### 22 nt
sam_in <- read_tsv(paste0(data_folder, "/", sample, "_q10_22.sam"), col_names = F)

out_sam <- sam_in %>% select(1:11) %>% group_by(X1) %>% mutate(MM = paste0("NH:i:", n()), n = n())
write_tsv(select(out_sam, -n), paste0(results_folder, "/",  sample, "_q10_22_filtered.sam"), , col_names=F)

out_sam_summary <- out_sam %>% ungroup() %>% distinct(X1, .keep_all = T) %>% group_by(n) %>% summarise(multi = n())
out_sam_summary

write_tsv(out_sam_summary, paste0(results_folder, "/", sample, "_q10_22_summary.tsv"))

##### 24 nt
sam_in <- read_tsv(paste0(data_folder, "/", sample, "_q10_24.sam"), col_names = F)

out_sam <- sam_in %>% select(1:11) %>% group_by(X1) %>% mutate(MM = paste0("NH:i:", n()), n = n())
write_tsv(select(out_sam, -n), paste0(results_folder, "/",  sample, "_q10_24_filtered.sam"), , col_names=F)

out_sam_summary <- out_sam %>% ungroup() %>% distinct(X1, .keep_all = T) %>% group_by(n) %>% summarise(multi = n())
out_sam_summary

write_tsv(out_sam_summary, paste0(results_folder, "/", sample, "_q10_24_summary.tsv"))
