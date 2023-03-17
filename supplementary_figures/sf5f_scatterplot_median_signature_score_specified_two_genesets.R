# Yige Wu @WashU Apr 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggpubr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input median signature scores per cluster
results_mat <- readRDS(file = "./Resources/Analysis_Results/signature_scores/process_signature_scores/make_median_signature_scores_per_res1_cluster/20220427.v1/Median_signature_scores_per_res1_cluster.20220427.v1.RDS")
## input the gene set auto-correlation results as a measure of consistency
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")

# make plot data -----------------------------------------------------------------
sigCorr_filtered_df <- sigCorr_df %>%
  filter(grepl(pattern = "HALLMARK", x = gene_set)) %>%
  filter(gene_set != "HALLMARK_UV_RESPONSE") %>%
  filter(FDR < 0.05 & C > 0.1) 
genesets_test <- sigCorr_filtered_df$gene_set
genesets_test <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
pairwise <- combn(genesets_test, 2)

pvalue_vec <- NULL
r_vec <- NULL
for (i in 1:ncol(pairwise)) {
  geneset_x <- pairwise[1,i]
  geneset_y <- pairwise[2,i]
  
  plotdata_df <- data.frame(t(results_mat[c(geneset_x, geneset_y),]))
  colnames(plotdata_df) <- c("x_plot", "y_plot")
  plotdata_df$cluster <- rownames(plotdata_df)
  
  # # plot --------------------------------------------------------------------
  # p <- ggscatter(data = plotdata_df, x = "x_plot", y = "y_plot",
  #                add = "reg.line",  # Add regressin line
  #                label = "cluster", font.label = c(14, "plain"),
  #                add.params = list(color = "grey", fill = "lightgray", linetype = 2)
  # )
  # p <- p + stat_cor(method = "pearson",
  #                   label.x = min(plotdata_df$x_plot),
  #                   label.y = (max(plotdata_df$y_plot)), size = 7)
  # # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " 1-month response"))
  # p <- p + theme_classic(base_size = 17)
  # p <- p + xlab(paste0(geneset_x, " score"))
  # p <- p + ylab(paste0(geneset_y, " score"))
  # # p <- p + xlim(c(min(plotdata_tmp_df$x_plot)-0.11, max(plotdata_tmp_df$x_plot)+0.11))
  # # p <- p + ylim(c(min(plotdata_tmp_df$y_plot)-0.05, max(plotdata_tmp_df$y_plot)+0.1))
  # p <- p + theme(axis.text = element_text(color = "black"))
  # file2write <- paste0(dir_out, geneset_x, ".vs.", geneset_y, ".png")
  # png(file2write, width = 600, height = 600, res = 150)
  # print(p)
  # dev.off()
  
  # cor_res <- cor(x = plotdata_df$x_plot, y = plotdata_df$y_plot, method = "pearson")
  # pvalue_vec <- c(pvalue_vec, cor_res)
  
  ## write source data
  write.table(x = plotdata_df, file = paste0("~/Desktop/SF5f.", geneset_x,".vs." , geneset_y, ".SourceData.tsv"), quote = F, sep = "\t", row.names = F)
}
# pearson_result_df <- data.frame(t(pairwise))
# colnames()


