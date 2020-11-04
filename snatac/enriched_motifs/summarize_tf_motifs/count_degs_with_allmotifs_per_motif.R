# Yige Wu @WashU Oct 2020

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
## input known TF relationship
tf2target_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/Transcriptional/omnipathdb.transcriptional.20200908.txt")
## input the DEG-TF matrix
deg2tf_wide_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/DEGs_with_TFs_inDARs.tsv")

# melt the DEG-TF table ---------------------------------------------------
deg2tf_long_df <- melt(deg2tf_wide_df)
deg2tf_long_filtered_df <- deg2tf_long_df %>%
  filter(value > 0) %>%
  rename(genesymbol_deg = V1) %>%
  rename(motifname = variable)

# count -------------------------------------------------------------------
count_degs_per_motif_df <- deg2tf_long_filtered_df %>%
  group_by(motifname) %>%
  summarise(count_degs = n()) %>%
  arrange(desc(count_degs))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Count_DEGs_with_Promoter_or_DistalMotif_per_Motif", ".tsv")
write.table(x = count_degs_per_motif_df, file = file2write, quote = F, sep = "\t", row.names = )



