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
  "ggplot2"
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
young_alltumors_df <- fread(data.table = F, input = "./Resources/Analysis_Results/process_external/process_Young_scRNA_Science_2018/findmarkers/findmarker_tumorcells_vs_nontumor_ccRCCpatients_tumorsonly/20220622.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
young_eachtumor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/process_external/process_Young_scRNA_Science_2018/findmarkers/findmarker_tumorcells_vs_nontumor_ccRCCpatients_byeachtumor/20220622.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
young_alltumors_vs_normal_df <- fread(data.table = F, input = "./Resources/Analysis_Results/process_external/process_Young_scRNA_Science_2018/findmarkers/findmarker_tumortissue_tumorcells_vs_normaltissue_nontumor_ccRCCpatients/20220622.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
bulk_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")

# specify parameters ---------------------------------------------------
genes_filter <- c("CA9", "SNAP25", "TGFA", "PLIN2", "ABCC3", "PHKA2", "KCTD3", "FTO", "SEMA6A", "EPHA6", "ABLM3", "PLEKHA1", "SLC6A3", "SHISA9", "CP", "PCSK6", "NDRG1", "EGFR", "ENPP3", "COL23A1", "UBE2D2")

# preprocess --------------------------------------------------------------
young_eachtumor_sum_df <- young_eachtumor_df %>%
  group_by(gene_symbol) %>%
  summarise(avg_log2FC = mean(avg_log2FC), p_val_adj = 10^(mean(log10(p_val_adj))))
genes_process_df <- merge(x = young_alltumors_df %>%
                            select(p_val_adj, avg_log2FC, gene_symbol), 
                          y = young_eachtumor_sum_df %>%
                            select(p_val_adj, avg_log2FC, gene_symbol),  by = c("gene_symbol"), all = T, suffixes = c(".young.alltumors", ".young.eachtumor"))
genes_process_df <- merge(x = genes_process_df,
                          y = young_alltumors_vs_normal_df %>%
                            rename(avg_log2FC.young.alltumors_vs_normal = avg_log2FC) %>%
                            select(gene_symbol, avg_log2FC.young.alltumors_vs_normal), by = c("gene_symbol"), all = T)
genes_process_df <- genes_process_df %>%
  mutate(ensembl_gene_id = paste0("ENSG", str_split_fixed(string = gene_symbol, pattern = "\\-ENSG", n = 2)[,2])) %>%
  mutate(gene_uniq_id = gene_symbol) %>%
  mutate(gene_symbol = str_split_fixed(string = gene_symbol, pattern = "\\-ENSG", n = 2)[,1])

genes_process_df <- merge(x = genes_process_df,
                          y = bulk_df %>%
                            rename(gene_symbol = Gene) %>%
                            rename(log2FC.wu_eachtumor = avg_log2FC.mean.TumorcellsvsNontumor) %>%
                            select(gene_symbol, log2FC.bulkRNA, log2FC.bulkRNA, log2FC.wu_eachtumor), by = c("gene_symbol"), all = T)

# make plot data ----------------------------------------------------------
plotdata_wide_df <- genes_process_df %>%
  filter(gene_symbol %in% genes_filter)
plotdata_df <- melt(plotdata_wide_df)
plotdata_df <- plotdata_df %>%
  mutate(data_type = ifelse(variable == "avg_log2FC.mean.TumorcellsvsNontumor", dataname_snrna,
                            ifelse(variable == "log2FC.snATAC", dataname_snatac, 
                                   ifelse(variable == "log2FC.bulkRNA", dataname_bulk_rna, dataname_bulk_protein)))) %>%
  mutate(foldchange = 2^value) %>%
  mutate(y_plot = Gene)
summary(plotdata_df$foldchange)
plotdata_df <- plotdata_df %>%
  mutate(x_plot = ifelse(foldchange >= 10, 10, foldchange))
plotdata_df$y_plot <- factor(x = plotdata_df$Gene, levels = plotdata_wide_df$Gene)
plotdata_df$data_type <- factor(x = plotdata_df$data_type, levels = c(dataname_snrna, dataname_snatac, dataname_bulk_rna, dataname_bulk_protein))

## make colors
display.brewer.all()
colors_datatype <- brewer.pal(n = 5, name = "Set1")[c(1, 3, 4, 5)]
names(colors_datatype) <- c(dataname_snrna, dataname_bulk_rna, dataname_bulk_protein, dataname_snatac)

# plot --------------------------------------------------------------------
p <- ggplot()
# p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = data_type, color = as.character(cell_group == "Macrophages")))
p <- p + geom_dotplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = data_type),
                      binaxis='y', stackdir='center', position=position_dodge(0.6), alpha = 0.7)
p <- p + scale_fill_manual(values = colors_datatype)
# p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
# p <- p + geom_hline(yintercept = 1, linetype = 2, alpha = 0.5)
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + scale_y_continuous(breaks = seq(0, 10, 2))
p <- p + ylab("Fold change")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
p <- p + theme(axis.text.y = element_text(size = 12, color = "black"), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12, color = "black"), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
p <- p + theme(legend.position = "top")
p <- p + guides(fill = guide_legend(override.aes = list(size=4), nrow = 4, title = NULL, label.theme = element_text(size = 12)))
file2write <- paste0(dir_out, "Foldchanges", ".png")
png(file2write, width = 600, height = 900, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "Foldchanges", ".pdf")
pdf(file2write, width = 4.25, height = 7, useDingbats = F)
print(p)
dev.off()



