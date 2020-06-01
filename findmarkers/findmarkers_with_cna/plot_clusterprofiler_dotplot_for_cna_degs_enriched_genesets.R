# Yige Wu @WashU May 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(clusterProfiler)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input ORA and GSEA results
oraresult_list <- readRDS(file = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/run_clusterprofiler_on_degs_associated_with_cna/20200519.v1/ORA.Result.20200519.v1.RDS")
gsearesult_list <- readRDS(file = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/run_clusterprofiler_on_degs_associated_with_cna/20200519.v1/GSEA.Result.20200519.v1.RDS")
## input known cnv info
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# plot dotplot for ORA result ---------------------------------------------
for (cna_gene_plot in names(oraresult_list)) {
  ## set output directory
  dir_out_sub <- paste0(dir_out, cna_gene_plot, "/")
  dir.create(dir_out_sub)
  ## get expected cnv states
  cna_gene_expected_state <- unique(knowncnvgenes_df$CNV_Type[knowncnvgenes_df$Gene_Symbol == cna_gene_plot])
  
  p <- dotplot(object = oraresult_list[[cna_gene_plot]])
  p <- p + ggtitle(label = paste0("ORA for Differentially Expressed Genes\nAssociated with ", cna_gene_plot, " ", cna_gene_expected_state))
  # p
  file2write <- paste0(dir_out_sub, cna_gene_plot, ".ORA.dotplot.png")
  png(filename = file2write, res = 150, width = 1500, height = 700)
  print(p)
  dev.off()
}

# plot dotplot for GEEA result ---------------------------------------------
for (cna_gene_plot in names(gsearesult_list)) {
  ## set output directory
  dir_out_sub <- paste0(dir_out, cna_gene_plot, "/")
  dir.create(dir_out_sub)
  ## get expected cnv states
  cna_gene_expected_state <- unique(knowncnvgenes_df$CNV_Type[knowncnvgenes_df$Gene_Symbol == cna_gene_plot])
  ## extract plot object
  object_plot <- gsearesult_list[[cna_gene_plot]]
  if (nrow(object_plot@result) > 0) {
    p <- dotplot(object = object_plot)
    p <- p + ggtitle(label = paste0("GSEA for Differentially Expressed Genes\nAssociated with ", cna_gene_plot, " ", cna_gene_expected_state))
    # p
    file2write <- paste0(dir_out_sub, cna_gene_plot, ".GSEA.dotplot.png")
    png(filename = file2write, res = 150, width = 1500, height = 700)
    print(p)
    dev.off()
  }
}
