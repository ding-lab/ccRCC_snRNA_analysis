# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/snATAC_Processed_Data/Gene_Activity/AverageGeneActivity_ByTumorPTLOHSubcluster.v2.20210917.tsv", data.table = F)
## input the genes to plot
genes4score_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_degs/filter_tumormanualcluster_EMT_degs/20210908.v1/EMTModuleDown.EpithelialPT_DEGs.Filtered.tsv")

# preprocess --------------------------------------------------------------
## identify clusters with sufficient cell number
genes4score <- genes4score_df$genesymbol_deg

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes4score) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "ATACGeneActivity\\.", replacement = "")) %>%
  dplyr::filter(!(grepl(x = id_bycluster_byaliquot, pattern = "LOH"))) %>%
  dplyr::filter(!(grepl(x = id_bycluster_byaliquot, pattern = "NA"))) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter(!(cluster_name %in% c("", "CNA")))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_byaliquot, value.var = "value")
plot_data_raw_mat <- as.matrix(plot_data_wide_df[,-1])
## add row names
rownames(plot_data_raw_mat) <- plot_data_wide_df$V1
## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- rownames(plot_data_raw_mat)
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# calculate geneset score -------------------------------------------------
score_vec <- colMeans(plot_data_mat[genes4score,])*100
colanno_df <- data.frame(cluster_name = colnames(plot_data_mat), score = score_vec)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "EpithelialScore.tsv")
write.table(x = colanno_df, file = file2write, quote = F, sep = "\t", row.names = F)
