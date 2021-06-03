# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(biomaRt)
library(karyoploteR)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
deg_5qgain_negbinom_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_negbinom_is5qgain_all_ccRCC_vs_pt_on_katmai/20210602.v1/negbinom.latent.Is.5q.Gain.logfc.threshold0.25.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
deg_negbinom_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_negbinom_all_ccRCC_vs_pt_on_katmai/20210602.v1/negbinom.logfc.threshold0.25.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")

# intersect ---------------------------------------------------------------
genes_up_5qgain_negbinom <- deg_5qgain_negbinom_df$genesymbol_deg[deg_5qgain_negbinom_df$p_val_adj < 0.05 & deg_5qgain_negbinom_df$avg_log2FC > 0]
genes_up_negbinom <- deg_negbinom_df$genesymbol_deg[deg_negbinom_df$p_val_adj < 0.05 & deg_negbinom_df$avg_log2FC > 0]
intersect(genes_up_5qgain_negbinom, genes_up_negbinom) %>% length()

genes_up_5qgain_negbinom <- deg_5qgain_negbinom_df$genesymbol_deg[deg_5qgain_negbinom_df$p_val_adj < 0.001 & deg_5qgain_negbinom_df$avg_log2FC > 0]
length(genes_up_5qgain_negbinom)
genes_up_negbinom <- deg_negbinom_df$genesymbol_deg[deg_negbinom_df$p_val_adj < 0.001 & deg_negbinom_df$avg_log2FC > 0]
length(genes_up_negbinom)
intersect(genes_up_5qgain_negbinom, genes_up_negbinom) %>% length()