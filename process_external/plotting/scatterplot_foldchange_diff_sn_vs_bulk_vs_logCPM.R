# Yige Wu @WashU Jun 2022

# set up libraries and output directory -----------------------------------
## set working directory
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
  "ggpubr",
  "ggrepel"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
young_t_tumorcells_vs_n_nontumor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/process_external/process_Young_scRNA_Science_2018/findmarkers/findmarker_tumortissue_tumorcells_vs_normaltissue_nontumor_ccRCCpatients/20220622.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
wu_t_tumorcells_vs_n_nontumor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findmarker_wilcox_tumortissue_tumorcells_vs_normaltissue_normalcells_34samples_katmai/20220622.v2/logfcthreshold.0.5.minpct.0.1.mindiffpct.0.1.tsv")
genes_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")
bulk.rna.df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_tumor_vs_NAT_deg_on_cptac_ccRCC_discovery_cases/20210419.v1/Tumor_vs_NAT.glmQLFTest.OutputTables.tsv")

# specify parameters ---------------------------------------------------
genes_filter <- c("CA9", genes_plot_df$Gene)
genes_highlight <- c("CA9", "SNAP25", "TGFA", "PLIN2", "ABCC3", "PHKA2", "KCTD3", "FTO", "SEMA6A", "EPHA6", "ABLM3", "PLEKHA1", "SLC6A3", "SHISA9", "CP", "PCSK6", "NDRG1", "EGFR", "ENPP3", "COL23A1", "UBE2D2")

# preprocess --------------------------------------------------------------
foldchanges_df <- merge(x = data.frame(gene_symbol = genes_filter),
                        y = young_t_tumorcells_vs_n_nontumor_df %>%
                          mutate(ensembl_gene_id = paste0("ENSG", str_split_fixed(string = gene_symbol, pattern = "\\-ENSG", n = 2)[,2])) %>%
                          mutate(gene_uniq_id = gene_symbol) %>%
                          mutate(gene_symbol = str_split_fixed(string = gene_symbol, pattern = "\\-ENSG", n = 2)[,1]) %>%
                          rename(log2FC.young = avg_log2FC) %>%
                          rename(fdr.young = p_val_adj) %>%
                          select(gene_symbol, log2FC.young, fdr.young), by = c("gene_symbol"), all.x = T)

foldchanges_df <- merge(x = foldchanges_df,
                        y = bulk.rna.df %>%
                          rename(gene_symbol = hgnc_symbol) %>%
                          rename(log2FC.bulkRNA = logFC) %>%
                          rename(fdr.bulkRNA = FDR) %>%
                          select(gene_symbol, log2FC.bulkRNA, logCPM, fdr.bulkRNA), by = c("gene_symbol"), all.x = T)
foldchanges_df <- merge(x = foldchanges_df,
                        y = wu_t_tumorcells_vs_n_nontumor_df %>%
                          rename(log2FC.wu = avg_log2FC) %>%
                          rename(fdr.wu = p_val_adj) %>%
                          select(gene_symbol, log2FC.wu, fdr.wu), by = c("gene_symbol"), all.x = T)

# log2FC.wu_vs_log2FC.bulkRNA --------------------------------------------------------------------
plotdata_df <- foldchanges_df %>%
  filter(!is.na(log2FC.wu) & !is.na(log2FC.bulkRNA) & !is.na(log2FC.young)) %>%
  filter(fdr.wu < 0.05) %>%
  mutate(x_plot = log2FC.bulkRNA) %>%
  mutate(y_plot = log2FC.wu) %>%
  mutate(diff.sn_vs_bulkrna = ifelse(abs(log2FC.wu - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.wu - log2FC.bulkRNA > 0 , "higher", "lower"))) %>%
  mutate(diff.sc_vs_bulkrna = ifelse(abs(log2FC.young - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.young - log2FC.bulkRNA > 0 , "higher", "lower")))
table(plotdata_df$diff.sn_vs_bulkrna)
116/243
27/243
100/243
table(plotdata_df$diff.sc_vs_bulkrna[plotdata_df$diff.sn_vs_bulkrna == "lower"])

p <- ggplot()
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey50")
p <- p + geom_abline(slope = 1, intercept = 1, linetype = 2, color = "grey50")
p <- p + geom_abline(slope = 1, intercept = -1, linetype = 2, color = "grey50")
p <- p + geom_point(data = subset(plotdata_df, diff.sn_vs_bulkrna == "higher"), mapping = aes(x = x_plot, y = y_plot), color = "red")
p <- p + geom_point(data = subset(plotdata_df, diff.sn_vs_bulkrna == "small diff"), mapping = aes(x = x_plot, y = y_plot), color = "black")
p <- p + geom_point(data = subset(plotdata_df, diff.sn_vs_bulkrna == "lower"), mapping = aes(x = x_plot, y = y_plot), color = "blue")
p <- p + geom_text_repel(data = subset(plotdata_df, diff.sn_vs_bulkrna == "lower" & gene_symbol %in% genes_highlight), 
                         mapping = aes(x = x_plot, y = y_plot, label = gene_symbol), max.overlaps = Inf)
p <- p + xlab(label = "Log2(bulk RNA-seq fold change)\nTumors vs. NATs")
p <- p + ylab(label = "Log2(snRNA-seq fold change)\nTumor cells vs. Non-tumor cells (NAT)")
p <- p + coord_fixed(ratio = 1)
p <- p + theme_classic(base_size = 12)
p
file2write <- paste0(dir_out, "log2FC.wu_vs_log2FC.bulkRNA", ".pdf")
pdf(file2write, width = 4, height = 4, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "log2FC.wu_vs_log2FC.bulkRNA", ".png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()

# log2FC.young_vs_log2FC.bulkRNA --------------------------------------------------------------------
plotdata_df <- foldchanges_df %>%
  filter(!is.na(log2FC.wu) & !is.na(log2FC.bulkRNA) & !is.na(log2FC.young)) %>%
  mutate(x_plot = log2FC.bulkRNA) %>%
  mutate(y_plot = ifelse(log2FC.young > 15, 15, ifelse(log2FC.young < -15, -15, log2FC.young))) %>%
  mutate(diff.sn_vs_bulkrna = ifelse(abs(log2FC.wu - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.wu - log2FC.bulkRNA > 0 , "higher", "lower"))) %>%
  mutate(diff.sc_vs_bulkrna = ifelse(abs(log2FC.young - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.young - log2FC.bulkRNA > 0 , "higher", "lower")))
table(plotdata_df$diff.sn_vs_bulkrna)
116/243
27/243
100/243
table(plotdata_df$diff.sc_vs_bulkrna[plotdata_df$diff.sn_vs_bulkrna == "lower"])

p <- ggplot()
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey50")
p <- p + geom_abline(slope = 1, intercept = 1, linetype = 2, color = "grey50")
p <- p + geom_abline(slope = 1, intercept = -1, linetype = 2, color = "grey50")
p <- p + geom_point(data = subset(plotdata_df, diff.sc_vs_bulkrna == "higher"), mapping = aes(x = x_plot, y = y_plot), color = "red")
p <- p + geom_point(data = subset(plotdata_df, diff.sc_vs_bulkrna == "small diff"), mapping = aes(x = x_plot, y = y_plot), color = "black")
p <- p + geom_point(data = subset(plotdata_df, diff.sc_vs_bulkrna == "lower"), mapping = aes(x = x_plot, y = y_plot), color = "blue")
p <- p + geom_text_repel(data = subset(plotdata_df, diff.sn_vs_bulkrna == "lower" & gene_symbol %in% genes_highlight), 
                         mapping = aes(x = x_plot, y = y_plot, label = gene_symbol), max.overlaps = Inf)
p <- p + xlab(label = "Log2(bulk RNA-seq fold change)\nTumors vs. NATs")
p <- p + ylab(label = "Log2(scRNA-seq fold change)\nTumor cells vs. Non-tumor cells (NAT)")
# p <- p + coord_fixed(ratio = 1)
p <- p + theme_classic(base_size = 12)
p
file2write <- paste0(dir_out, "log2FC.young_vs_log2FC.bulkRNA", ".pdf")
pdf(file2write, width = 4, height = 4, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "log2FC.young_vs_log2FC.bulkRNA", ".png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()

# log2FC.wu_-_log2FC.bulkRNA --------------------------------------------------------------------
plotdata_df <- foldchanges_df %>%
  filter(!is.na(log2FC.wu) & !is.na(log2FC.bulkRNA)) %>%
  mutate(x_plot = logCPM) %>%
  mutate(y_plot = log2FC.wu - log2FC.bulkRNA) %>%
  mutate(diff.sn_vs_bulkrna = ifelse(abs(log2FC.wu - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.wu - log2FC.bulkRNA > 0 , "higher", "lower"))) %>%
  mutate(diff.sc_vs_bulkrna = ifelse(abs(log2FC.young - log2FC.bulkRNA) < 1 , "small diff",
                                     ifelse(log2FC.young - log2FC.bulkRNA > 0 , "higher", "lower")))

p <- ggscatter(data = plotdata_df, x = "x_plot", y = "y_plot",
               add = "reg.line",
               add.params = list(color = "grey", fill = "lightgray", linetype = 2)
)
p <- p + stat_cor(method = "pearson",
                  label.x = min(plotdata_df$x_plot),
                  label.y = (max(plotdata_df$y_plot)), size = 7)
p <- p + geom_text_repel(data = subset(plotdata_df, diff.sn_vs_bulkrna == "lower" & gene_symbol %in% genes_highlight), 
                         mapping = aes(x = x_plot, y = y_plot, label = gene_symbol), max.overlaps = Inf, min.segment.length = 0, force = 2)
p <- p + ylab(label = "log2FC (snRNA-seq) - log2FC (bulk RNA-seq)")
p <- p + xlab(label = "log2CPM")

file2write <- paste0(dir_out, "diff.sn_bulkrna_vs_CPM", ".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "diff.sn_bulkrna_vs_CPM", ".pdf")
pdf(file2write, width = 4.25, height = 7, useDingbats = F)
print(p)
dev.off()


