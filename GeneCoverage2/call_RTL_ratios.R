call_RTL_ratios <- function(Sample, dataFolder, outDir, readlength, min_intergenic_distance){

  # read length should be doubled if PE eg 2 x 100 = 200 bp

  ####### args for debugging
# Sample = "Sample_alx8_277_9"
# dataFolder = "PerGeneCoverageBinned/per_gene_tables"
# outDir = "PerGeneCoverageBinned/RTL_ratios"
# readlength = 200
# min_intergenic_distance = 50

  ####### set up

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

####### process
print(paste0("Processing ", Sample))
# read in data
per_gene_bin_cov <- read_csv(paste0(dataFolder, "/", Sample, "_gene_cov_bin100.csv"))

total_input_genes <- nrow(distinct(per_gene_bin_cov, gene))

print(paste0("total input genes ", total_input_genes))

########## 1. Remove genes < 50 bp intergenic data

### testing
# per_gene_bin_cov_bk <- per_gene_bin_cov
# per_gene_bin_cov <- slice(per_gene_bin_cov_bk, 1:100000)

# get list of genes with > 50 bp intergenic
# if there are 5 bins then there should be 50bp of intergenic distance because distance file had 0 if no data and all position in genome were assigned
# note that the input file I am using I ran bedtools closest but really i could have just found the nearest upstream gene and not worry about the downstream...
min_number_bins <- min_intergenic_distance/10
# min_number_bins <- 3
print(paste0("min downstream intergenic distance ", min_intergenic_distance, "; min number of bins ", min_number_bins))

genes_50bp_intergenic <- per_gene_bin_cov %>%
  group_by(gene, gene_cat, side) %>%
  summarise(n = n()) %>%
  filter(gene_cat == "down" & side == "Sense" & n >= min_number_bins) %>%
  ungroup() %>%
  select(gene)

print(paste0("genes passing ", min_intergenic_distance, " bp intergenic filter ", nrow(genes_50bp_intergenic)))

######### 2. [Remove genes < 0.03 CPM]

# calculate CPM ish
 total_counts <- per_gene_bin_cov %>%
  filter(gene_cat == "genic" & side == "Sense") %>%
  pull(sum_bin_coverage)

 # scaling factor is to estimate CPM therefore / total counts by 1M * read length to get an estimate
  scaling_factor <- sum(as.numeric(total_counts), na.rm = T) /(1000000*readlength)
scaling_factor

 genes_CPM <- per_gene_bin_cov %>%
  filter(gene_cat == "genic" & side == "Sense") %>%
   group_by(gene) %>%
   summarise(total_gene_counts = sum(as.numeric(sum_bin_coverage), na.rm = T)) %>%
  mutate(psudo_CPM = total_gene_counts/scaling_factor) %>%
   ungroup() %>%
   select(gene, psudo_CPM)

 # print(paste0("writing psudo_CPM file"))
 #
 # write_csv(genes_CPM, paste0(outDir, "/", Sample, "_psudo_CPM.csv"))

 # combine filter lists
 # genes_passing_filtering <- inner_join(genes_50bp_intergenic, genes_0.3CPM, by = "gene")
 #
 # print(paste0("genes passing 0.03 CPM filter and 50 bp intergenic filter ", nrow(genes_passing_filtering)))

######### Filter, Subset to 500 bp +/- TSS

per_gene_bin_cov_filter <- per_gene_bin_cov %>%
  filter(gene %in% genes_50bp_intergenic$gene & side == "Sense" & gene_cat %in% c("genic", "down")) %>%
  left_join(bin_key, by = c("gene_cat", "bin100")) %>%
  arrange(side, bin_name_fct) %>%
  # distinct(bin100) %>%
  # pull(bin100)
  filter(gene_cat == "genic" & bin100 >= 51 | gene_cat == "down" & bin100 <= 50) %>%
  group_by(gene, gene_cat) %>%
  summarise(mean_cov = mean(sum_bin_coverage, na.rm = T)) %>%
  spread(key = gene_cat, value = mean_cov) %>%
  # filter(down >= 2) %>% # Remove genes <2x mean intergenic coverage
  mutate(RT_ratio = log2((down+0.25)/(genic+0.25))) %>% # include variance norm to account for zeros
  # filter(RT_ratio < log2(8)) # remove mean intergenic coverage >8 times than of genic coverage
  left_join(genes_CPM, by = "gene")

per_gene_bin_cov_filter

per_gene_bin_cov_filter %>% filter(!is.finite(RT_ratio))

 print(paste0("total genes with calculated read-through ratios ", nrow(per_gene_bin_cov_filter)))

 write_csv(per_gene_bin_cov_filter, paste0(outDir, "/", Sample, "_RTL_ratios.csv"))
}
