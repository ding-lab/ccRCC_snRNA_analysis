# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210610.v1/BAP1_DEGs.United.snRNA.bulkRNA.Protein.20210610.v1.tsv")

# plot by cutoff ----------------------------------------------------------
plotdata_df <- NULL
for (cutoff_tmp in 1) {
# for (cutoff_tmp in 1:5) {
  plotdata_tmp_df <- degs_df %>%
    filter(FDR.cnvcorrected < 0.05 & !is.na(FDR.cnvcorrected)) %>%
    filter(!is.na(Num_up)) %>%
    mutate(ifelse())
    filter((Num_up == 0 & Num_sig_down >= cutoff_tmp) | (Num_down == 0 & Num_sig_up >= cutoff_tmp)) %>%
    select(BAP1_vs_OtherTumor_snRNA) %>%
    table() %>%
    as.data.frame() %>%
    rename(BAP1_vs_OtherTumor_snRNA = '.') %>%
    mutate(min_sig_samples = cutoff_tmp)
  plotdata_df <- rbind(plotdata_df, plotdata_tmp_df)
}
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = min_sig_samples, y = Freq, fill = BAP1_vs_OtherTumor_snRNA), stat = "identity")
p