####### metaplot function

# Gene_list_meta_plot_PerGeneCoverageBinned <- function(Sample, dataFolder, filter_list_name, gene_list_path){

  # This function requires thet the gene list has a column called "Gene_ID" with TAIR ATGs with no trailing "." isoform numbers
  # This function should be run from the folder containing the folder "PerGeneCoverageBinned"

  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  # trailingOnly=TRUE means that only your arguments are returned
  Sample <- args[1]
  Sample
  dataFolder <- args[2]
  dataFolder
  filter_list_name <- args[3]
  filter_list_name
  gene_list_path <- args[4]
  gene_list_path
  library_layout <- args[5]
  library_layout

library(tidyverse)

###### Args
# Sample = "Sample_alx8_277_9"
# dataFolder = "PerGeneCoverageBinned/per_gene_tables"
# filter_list_name <- "RTLs"
# gene_list_path <- "~/ws/refseqs/TAIR10/RTLs.csv"

##### set up
outDir = paste0(dataFolder, "_", filter_list_name)
dir.create(outDir, showWarnings = F)

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

###### read in data
per_gene_bin_cov <- read_csv(paste0(dataFolder, "/", Sample, "_gene_cov_bin100.csv"))

per_gene_bin_cov <- per_gene_bin_cov %>%
  separate(gene, into = c("Gene_ID", "variant"))

####### gene list filter
# read in gene list
# gene_list <- read_csv("RTLs.csv", skip = 4)
gene_list <- read_csv(gene_list_path)

colnames(gene_list)

per_gene_bin_cov_filtered <- per_gene_bin_cov %>%
  filter(Gene_ID %in% gene_list$Gene_ID)

####### summarise for metaplot

bin_cov_summary <- per_gene_bin_cov_filtered %>%
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


write_csv(bin_cov_summary, paste0(outDir, "/", Sample, "_gene_cov_bin100_meta.csv"))

# metaplot sum
pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = meta_sum_bin_coverage, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

# metaplot sum then log
pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot_log.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = log_cov, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

# metaplot mean(mean(log))
pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot_mean_log.pdf"))
g <- ggplot(bin_cov_summary, aes(x = position, y = meta_mean_bin_coverage, colour = side)) +
  geom_line() +
  theme_classic()
print(g)
dev.off()

}else{
  if(library_layout == "nonstranded"){
  print("library_layout = nonstranded")

  ###### read in data
  per_gene_bin_cov <- read_csv(paste0(dataFolder, "/", Sample, "_gene_cov_bin100.csv"))

  per_gene_bin_cov <- per_gene_bin_cov %>%
    separate(gene, into = c("Gene_ID", "variant"))

  ####### gene list filter
  # read in gene list
  # gene_list <- read_csv("RTLs.csv", skip = 4)
  gene_list <- read_csv(gene_list_path)

  colnames(gene_list)

  per_gene_bin_cov_filtered <- per_gene_bin_cov %>%
    filter(Gene_ID %in% gene_list$Gene_ID)

  ####### summarise for metaplot

  bin_cov_summary <- per_gene_bin_cov_filtered %>%
    group_by(gene_cat, bin100) %>%
    summarise(meta_sum_bin_coverage = sum(sum_bin_coverage),
    meta_mean_bin_coverage = mean(mean_bin_coverage)) %>%
    ungroup() %>%
    left_join(bin_key, by = c("gene_cat", "bin100")) %>%
    arrange(bin_name_fct) %>%
    mutate(position = row_number(),
           log_cov = ifelse(meta_sum_bin_coverage == 0, 0,
                            ifelse(meta_sum_bin_coverage > 0, log(meta_sum_bin_coverage),
                            ifelse(meta_sum_bin_coverage < 0, (log(abs(meta_sum_bin_coverage)))*-1, 0))))


  write_csv(bin_cov_summary, paste0(outDir, "/", Sample, "_gene_cov_bin100_meta.csv"))

  # metaplot sum
  pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot.pdf"))
  g <- ggplot(bin_cov_summary, aes(x = position, y = meta_sum_bin_coverage)) +
    geom_line() +
    theme_classic()
  print(g)
  dev.off()

  # metaplot sum then log
  pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot_log.pdf"))
  g <- ggplot(bin_cov_summary, aes(x = position, y = log_cov)) +
    geom_line() +
    theme_classic()
  print(g)
  dev.off()

  # metaplot mean(mean(log))
  pdf(paste0(outDir, "/", Sample, "_gene_cov_bin100_metaplot_mean_log.pdf"))
  g <- ggplot(bin_cov_summary, aes(x = position, y = meta_mean_bin_coverage)) +
    geom_line() +
    theme_classic()
  print(g)
  dev.off()

  }else{
  print("strandedness not specified")
  }
  }

#}
