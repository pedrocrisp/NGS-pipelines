# # Created by Peter Crisp
# March 2018
# pedrocrisp at gmail.com

# Gene coverage bining.

# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
library(tidyverse)

# outDir
outDir = "PerGeneCoverageBinned"

args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
Sample <- args[1]
Sample
beds_folder <- args[2]
beds_folder
library_layout <- args[3]
library_layout

##########
# debug
# Sample = "test"
# beds_folder = "test_stranded_cov_beds"
# # Sample = "Salt_heat_2"
# # beds_folder = "tdf_for_igv_coverage_beds_nonstranded"
# outDir = "test_GeneCoverage_bined"

# Sample = "Sample_alx8_277_9"
# outDir = "GeneCoverage_bined"
# beds_folder = "tdf_for_igv_coverage_beds_full"

print(Sample)

sPath <- paste0(beds_folder, "/", Sample, "/")

#meta plots out
dir.create(outDir, showWarnings = F, recursive = T)

outFolder_plots <- paste0(outDir, "/meta_plots")
dir.create(outFolder_plots, showWarnings = F, recursive = T)
#meta tables out
outFolder_meta_tables <- paste0(outDir, "/meta_table")
dir.create(outFolder_meta_tables, showWarnings = F, recursive = T)
# per_gene_tables out
outFolder_per_gene_tables <- paste0(outDir, "/per_gene_tables")
dir.create(outFolder_per_gene_tables, showWarnings = F, recursive = T)

##########

print("create breaks key")
 # get breaks for up and down stream this will mean in even breaks and appropriate missing values in metaplot or histogram
 # note: at a later point I could use the mids to label a plot properly
 mids_up <- hist(c(-1:-1000), breaks = 100, plot = F)[["mids"]]
 breaks_up <- hist(c(-1:-1000), breaks = 100, plot = F)[["breaks"]]
 mids_down <- hist(c(1:1000), breaks = 100, plot = F)[["mids"]]
 breaks_down <- hist(c(1:1000), breaks = 100, plot = F)[["breaks"]]

 mids_gene <- hist(c(0:100), breaks = 100, plot = F)[["mids"]]

 bin_key_up <- tibble(bin100 = c(1:100), gene_cat = "up", bin_name = mids_up)
 bin_key_gene <- tibble(bin100 = c(1:100), gene_cat = "genic", bin_name = mids_gene)
 bin_key_down <- tibble(bin100 = c(1:100), gene_cat = "down", bin_name = mids_down)

bin_key <- bind_rows(bin_key_up, bin_key_gene, bin_key_down) %>%
  mutate(bin_name_fct = paste0(bin_name, "_", gene_cat)) %>%
  mutate(bin_name_fct = factor(bin_name_fct, levels = bin_name_fct))

if(library_layout == "stranded"){
  print("library_layout = stranded")

############################## plus bed
print("process plus bed")
# read in required columns
coverage_file <- read_tsv(paste0(beds_folder, "/", Sample, "/", Sample, ".plus.dist.1k.bed"),
                          col_names = c("Chr",	"start",	"stop",	"coverage",	"feature_Chr",	"feature_start",	"feature_stop",	"gene",	"score",	"strand",	"annotation",	"feature",	"feature_number",	"gene_start",	"gene_stop",	"Gene_model_number",	"Distance"),
                          cols_only(
                              Chr = col_character(),
                              coverage = col_double(),
                              gene = col_character(),
                              strand = col_character(),
                              feature = col_character(),
                              Distance = col_integer()
                            ))

# get IDs for genes shorter than 300 bp of exon
# there are none?
print("get IDs for genes shorter than 300 bp of exon")
 short_genes <- coverage_file %>%
  filter(feature == "exon" | Distance == 0) %>%
  group_by(gene) %>%
  filter(n() < 300) %>%
   distinct(gene) %>%
   pull(gene)

########## split + and -
print("split + and -, remove short genes and organelles")
coverage_file_plus <- coverage_file %>% filter(strand == "+") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)
coverage_file_minus <- coverage_file %>% filter(strand == "-") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)

######## bin genic region and upstream and downstream

####### + strand genes

 # bin data per locus
print("bin data per locus + strand")
per_gene_bin_cov_plus <- coverage_file_plus %>%
  filter(feature == "exon") %>% # remove introns
  mutate(real_dist = -1*Distance, side = "Sense") %>%
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!gene %in% short_genes) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F, include.lowest = T),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
  mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100, side) %>%
  # filter(gene_cat =="genic")
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum and mean coverage for each bin
  ungroup()

per_gene_bin_cov_plus %>% filter(sum_bin_coverage > 0)
per_gene_bin_cov_plus %>% filter(sum_bin_coverage < 0)
per_gene_bin_cov_plus %>% filter(sum_bin_coverage == 0)
per_gene_bin_cov_plus

# bin_cov_summary <- per_gene_bin_cov_plus %>%
#   group_by(gene_cat, bin100, side) %>%
#   summarise(meta_sum_bin_coverage = sum(sum_bin_coverage)) %>%
#   ungroup() %>%
#   left_join(bin_key, by = c("gene_cat", "bin100")) %>%
#   arrange(side, bin_name_fct) %>%
#   group_by(side) %>%
#   mutate(position = row_number(),
#          log_cov = ifelse(meta_sum_bin_coverage == 0, 0,
#                           ifelse(meta_sum_bin_coverage > 0, log(meta_sum_bin_coverage),
#                                  ifelse(meta_sum_bin_coverage < 0, (log(abs(meta_sum_bin_coverage)))*-1, 0))))
#
# g <- ggplot(bin_cov_summary, aes(x = position, y = meta_sum_bin_coverage)) +
#   geom_line() +
#   theme_classic()
# print(g)

####### - strand genes

 # bin data per locus
print("bin data per locus - strand")
per_gene_bin_cov_minus <- coverage_file_minus %>%
  filter(feature == "exon") %>% # remove introns
  arrange(-row_number()) %>%
  mutate(real_dist = Distance, side = "Antisense", coverage = coverage *-1) %>% # Correct bedtools distance for strand of gene
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!gene %in% short_genes) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
  mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100, side) %>%
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum coverage for each bin
  ungroup()

per_gene_bin_cov_minus %>% filter(sum_bin_coverage > 0)
per_gene_bin_cov_minus %>% filter(sum_bin_coverage < 0)

####### combine + and -
per_gene_bin_cov_plus_bed <- bind_rows(per_gene_bin_cov_plus, per_gene_bin_cov_minus)

##############

# bin_cov_summary <- per_gene_bin_cov_plus_bed %>%
#   group_by(gene_cat, bin100, side) %>%
#   summarise(meta_sum_bin_coverage = sum(sum_bin_coverage)) %>%
#   ungroup() %>%
#   left_join(bin_key, by = c("gene_cat", "bin100")) %>%
#   arrange(side, bin_name_fct) %>%
#   group_by(side) %>%
#   mutate(position = row_number(),
#          log_cov = ifelse(meta_sum_bin_coverage == 0, 0,
#                           ifelse(meta_sum_bin_coverage > 0, log(meta_sum_bin_coverage),
#                                  ifelse(meta_sum_bin_coverage < 0, (log(abs(meta_sum_bin_coverage)))*-1, 0))))
#
# g <- ggplot(bin_cov_summary, aes(x = position, y = meta_sum_bin_coverage, colour = side)) +
#   geom_line() +
#   theme_classic()
# print(g)


############################## minus bed
print("process minus bed")
# read in required columns
coverage_file <- read_tsv(paste0(beds_folder, "/", Sample, "/", Sample, ".minus.dist.1k.bed"),
                          col_names = c("Chr",	"start",	"stop",	"coverage",	"feature_Chr",	"feature_start",	"feature_stop",	"gene",	"score",	"strand",	"annotation",	"feature",	"feature_number",	"gene_start",	"gene_stop",	"Gene_model_number",	"Distance"),
                          cols_only(
                              Chr = col_character(),
                              coverage = col_double(),
                              gene = col_character(),
                              strand = col_character(),
                              feature = col_character(),
                              Distance = col_integer()
                            ))

# get IDs for genes shorter than 300 bp of exon
# there are none?
print("get IDs for genes shorter than 300 bp of exon")
 short_genes <- coverage_file %>%
  filter(feature == "exon" | Distance == 0) %>%
  group_by(gene) %>%
  filter(n() < 300) %>%
   distinct(gene) %>%
   pull(gene)

########## split + and -
print("split + and -, remove short genes and organelles")
coverage_file_plus <- coverage_file %>% filter(strand == "+") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)
coverage_file_minus <- coverage_file %>% filter(strand == "-") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)

######## bin genic region and upstream and downstream

####### + strand genes

 # bin data per locus
print("bin data per locus + strand")
per_gene_bin_cov_plus <- coverage_file_plus %>%
  filter(feature == "exon") %>% # remove introns
  mutate(real_dist = -1*Distance, side = "Antisense") %>%
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!gene %in% short_genes) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100, side) %>%
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum coverage for each bin
  ungroup()

per_gene_bin_cov_plus %>% filter(sum_bin_coverage < 0)
per_gene_bin_cov_plus %>% filter(sum_bin_coverage > 0)

####### - strand genes

 # bin data per locus
print("bin data per locus - strand")
per_gene_bin_cov_minus <- coverage_file_minus %>%
  filter(feature == "exon") %>% # remove introns
  arrange(-row_number()) %>%
  mutate(real_dist = Distance, side = "Sense", coverage = coverage *-1) %>% # Correct bedtools distance for strand of gene
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!gene %in% short_genes) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
  mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100, side) %>%
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum coverage for each bin
  ungroup()

per_gene_bin_cov_minus %>% filter(sum_bin_coverage > 0)

####### combine + and -
per_gene_bin_cov_minus_bed <- bind_rows(per_gene_bin_cov_plus, per_gene_bin_cov_minus)



####### combine plus and minus beds
print("recombine summarise and write files")
per_gene_bin_cov <- bind_rows(per_gene_bin_cov_plus_bed, per_gene_bin_cov_minus_bed)

write_csv(per_gene_bin_cov, paste0(outFolder_per_gene_tables, "/", Sample, "_gene_cov_bin100.csv"))



####### summarise for metaplot

bin_cov_summary <- per_gene_bin_cov %>%
  group_by(gene_cat, bin100, side) %>%
  summarise(meta_sum_bin_coverage = sum(sum_bin_coverage),
  meta_mean_bin_coverage = mean(mean_bin_coverage)) %>%
  ungroup() %>%
  left_join(bin_key, by = c("gene_cat", "bin100")) %>%
  arrange(side, bin_name_fct) %>%
  group_by(side) %>%
  mutate(position = row_number(),
         log_cov = ifelse(meta_sum_bin_coverage == 0, 0,
                          ifelse(meta_sum_bin_coverage > 0, log(meta_sum_bin_coverage),
                          ifelse(meta_sum_bin_coverage < 0, (log(abs(meta_sum_bin_coverage)))*-1, 0))))


write_csv(bin_cov_summary, paste0(outFolder_meta_tables, "/", Sample, "_gene_cov_bin100_meta.csv"))

# metaplot sum
pdf(paste0(outFolder_plots, "/", Sample, "_gene_cov_bin100_metaplot.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = meta_sum_bin_coverage, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

# metaplot sum then log
pdf(paste0(outFolder_plots, "/", Sample, "_gene_cov_bin100_metaplot_log.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = log_cov, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

# metaplot mean(mean(log))
pdf(paste0(outFolder_plots, "/", Sample, "_gene_cov_bin100_metaplot_mean_log.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = meta_mean_bin_coverage, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

}else{
  if(library_layout == "nonstranded"){
  print("library_layout = nonstranded")

  #UPDATE TO INCLUDE NEW CODE FOR MEAN LOG IN STRANDED CODE ABOVE!!!
  # This has been added but is untested - use with caution and make sure to QC

  # read in required columns
coverage_file <- read_tsv(paste0(beds_folder, "/", Sample, "/", Sample, ".dist.1k.bed"),
                          col_names = c("Chr",	"start",	"stop",	"coverage",	"feature_Chr",	"feature_start",	"feature_stop",	"gene",	"score",	"strand",	"annotation",	"feature",	"feature_number",	"gene_start",	"gene_stop",	"Gene_model_number",	"Distance"),
                          cols_only(
                              Chr = col_character(),
                              coverage = col_double(),
                              gene = col_character(),
                              strand = col_character(),
                              feature = col_character(),
                              Distance = col_integer()
                            ))

# get IDs for genes shorter than 300 bp of exon
# there are none?
print("get IDs for genes shorter than 300 bp of exon")
 short_genes <- coverage_file %>%
  filter(feature == "exon" | Distance == 0) %>%
  group_by(gene) %>%
  filter(n() < 300) %>%
   distinct(gene) %>%
   pull(gene)

########## split + and -
print("split + and -, remove short genes and organelles")
coverage_file_plus <- coverage_file %>% filter(strand == "+") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)
coverage_file_minus <- coverage_file %>% filter(strand == "-") %>%
  filter(!gene %in% short_genes) %>% #remove short gene
  filter(!Chr %in% c("ChrC", "ChrM")) %>%
  select(-Chr)

######## bin genic region and upstream and downstream

####### plus

 # bin data per locus
print("bin data per locus plus strand")
per_gene_bin_cov_plus <- coverage_file_plus %>%
  filter(feature == "exon") %>% # remove introns
  mutate(real_dist = Distance) %>%
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!n() < 100) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
  mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100) %>%
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum coverage for each bin
  ungroup()

####### minus

 # bin data per locus
print("bin data per locus minus strand")
per_gene_bin_cov_minus <- coverage_file_minus %>%
  filter(feature == "exon") %>% # remove introns
  arrange(-row_number()) %>%
  mutate(real_dist = ifelse(strand=='+',-1*Distance,Distance)) %>% # Correct bedtools distance for strand of gene
  mutate(gene_cat = ifelse(real_dist < 0, "up", ifelse(real_dist == 0, "genic", "down"))) %>%
  mutate(gene_cat = factor(gene_cat, levels = c("up", "genic", "down"))) %>%
  group_by(gene, gene_cat) %>% # group
  filter(!n() < 100) %>% # remove any gene with total exon length < 300
  mutate(exon_position = row_number()) %>% # add numerical order number
  mutate(bin100 =
           ifelse(gene_cat == "genic", cut(exon_position, breaks = 100, labels = F),
           ifelse(gene_cat == "up", cut(real_dist, breaks = breaks_up, labels = F),
           ifelse(gene_cat == "down", cut(real_dist, breaks = breaks_down, labels = F), 0)))) %>% # bin the data
  mutate(log_cov = ifelse(coverage == 0, coverage, log(abs(coverage)))) %>%
  group_by(gene, gene_cat, bin100) %>%
  summarise(sum_bin_coverage = sum(coverage),
  mean_bin_coverage = mean(log_cov, na.rm = T)) %>% # sum coverage for each bin
  ungroup()

####### combine + and -
print("recombine summarise and write files")
per_gene_bin_cov <- bind_rows(per_gene_bin_cov_plus, per_gene_bin_cov_minus)

write_csv(per_gene_bin_cov, paste0(outFolder_gene_tables, "/", Sample, "_gene_cov_bin100.csv"))



####### summarise for metaplot

bin_cov_summary <- per_gene_bin_cov %>%
  group_by(gene_cat, bin100) %>%
  summarise(meta_sum_bin_coverage = sum(sum_bin_coverage),
  meta_mean_bin_coverage = mean(mean_bin_coverage)) %>%
  ungroup() %>%
  mutate(position = row_number())

write_csv(bin_cov_summary, paste0(outFolder_meta_tables, "/", Sample, "_gene_cov_bin100_meta.csv"))

# metaplot
pdf(paste0(outFolder_plots, "/", Sample, "_gene_cov_bin100_metaplot.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = log10(meta_sum_bin_coverage))) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

# metaplot mean(mean(log))
pdf(paste0(outFolder_plots, "/", Sample, "_gene_cov_bin100_metaplot_mean_log.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = meta_mean_bin_coverage, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

  }else{
  print("strandedness not specified")
}
}
