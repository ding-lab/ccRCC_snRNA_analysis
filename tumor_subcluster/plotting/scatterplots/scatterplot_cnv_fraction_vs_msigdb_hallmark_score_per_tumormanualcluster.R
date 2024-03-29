# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ggpubr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
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
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20201207.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20201207.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input gene sets to test
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")

# set parameters ----------------------------------------------------------
## specify the CNV gene and type to plot
cnvs_plot_df <- knowncnvgenes_df
## specify the pathway scores to plot
genesetnames_plot <- genesets_test_df$Description

# make data frame for plotting --------------------------------------------
## preprocess
scores_df <- scores_df %>%
  mutate(tumor_subcluster = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
dir_out1 <- paste0(dir_out, "pdf/"); dir.create(dir_out1)
dir_out1 <- paste0(dir_out, "png/"); dir.create(dir_out1)

for (genesetname in genesetnames_plot) {
  ## extract pathway score
  geneset_colname <- paste0(gsub(x = genesetname, pattern = "HALLMARK_", replacement = ""), "_Score")
  scores_tmp_df <- scores_df[, c("tumor_subcluster", geneset_colname)] 
  colnames(scores_tmp_df) <- c("tumor_subcluster", "geneset_score")
  for (i in 1:nrow(cnvs_plot_df)) {
    ## extract CNV data
    gene_cnv <- cnvs_plot_df$Gene_Symbol[i]
    cnv_type <- cnvs_plot_df$CNV_Type[i]
    dir_out2 <- paste0(dir_out1, gene_cnv, "_", cnv_type, "/"); dir.create(dir_out2)
    
    cnv_tmp_df <- cnv_3state_count_aliquots %>%
      filter(gene_symbol == gene_cnv) %>%
      filter(cna_3state == cnv_type)
    if (nrow(cnv_tmp_df) == 0) {
      next()
    }
    file2write <- paste0(dir_out2, gene_cnv, "_", cnv_type, ".", genesetname, ".pdf")
    # if (file.exists(file2write)) {
    #   next()
    # }
    ## make data frame for the NA data
    cluster_nonna_df <- cnv_3state_count_aliquots %>%
      filter(gene_symbol %in% gene_cnv) %>%
      select(tumor_subcluster) %>%
      unique() %>%
      mutate(Data_detected = T)
    ## merge
    plotdata_df <- merge(x = cnv_tmp_df, y = scores_tmp_df, by = c("tumor_subcluster"), all.y = T)
    plotdata_df <- merge(x = plotdata_df, y = cluster_nonna_df, by = c("tumor_subcluster"), all.y = T)
    ## filter
    plotdata_df <- plotdata_df %>%
      filter(!is.na(Data_detected) & Data_detected) %>%
      mutate(x_plot = ifelse(is.na(Fraction), 0, Fraction)) %>%
      mutate(y_plot = geneset_score) %>%
      mutate(sample_id = str_split_fixed(string = tumor_subcluster, pattern = "_", n = 2)[,1]) %>%
      filter(!is.na(x_plot) & !is.na(y_plot))
     
    if (nrow(plotdata_df) >= 5) {
      ## plot
      p <- ggscatter(data = plotdata_df, x = "x_plot", y = "y_plot",
                     add = "reg.line",  # Add regressin line
                     # label = "cluster", font.label = c(14, "plain"),
                     add.params = list(color = "grey", fill = "lightgray", linetype = 2)
      )
      p <- p + stat_cor(method = "pearson",
                        label.x = min(plotdata_df$x_plot),
                        label.y = quantile(plotdata_df$y_plot, 0.9), size = 5)
      p <- p + theme_classic()
      p <- p + ylab(label = geneset_colname) + xlab(paste0("% ", gene_cnv, " ", cnv_type, "\nper tumor cluster"))
      p <- p + theme(axis.text = element_text(color = "black"))
      file2write <- paste0(dir_out2, gene_cnv, "_", cnv_type, ".", genesetname, ".png")
      png(file2write, width = 600, height = 600, res = 150)
      print(p)
      dev.off()
      # file2write <- paste0(dir_out1, gene_cnv, "_", cnv_type, ".", genesetname, ".pdf")
      # pdf(file2write, width = 2.5, height = 2.5, useDingbats = F)
      # print(p)
      # dev.off()
    }
  }
}




