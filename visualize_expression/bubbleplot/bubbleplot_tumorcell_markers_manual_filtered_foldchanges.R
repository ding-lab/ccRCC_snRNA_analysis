# Yige Wu @WashU Sep 2021

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
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")

# specify parameters ---------------------------------------------------
## process genes to plot 
genes_filter <- genes_process_df$Gene
genes_filter <- genes_filter[!(genes_filter %in% c("PIK3CB", "ARHGEF28", "PTGER3", "PARD3", "GNG12", "EFNA5", "SPIRE1", "LIFR", "PKP4", "SORBS1", "PTPRM", "FBXO16", "PAM"))]
genes_filter <- genes_filter[!(genes_filter %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM"))]
## make 
dataname_snrna <- "Tumor cells vs. non-tumor cells\n(snRNA-seq)"
dataname_snrna <- "Tumor cells vs. non-tumor cells (snRNA-seq)"
dataname_bulk_rna <- "Tumors vs. NATs (bulk RNA-seq)"
dataname_bulk_protein <- "Tumors vs. NATs (bulk proteomics)"

# make plot data ----------------------------------------------------------
plotdata_wide_df <- genes_process_df %>%
  filter(Gene %in% genes_filter) %>%
  select(Gene, avg_log2FC.allTumorcellsvsPT, log2FC.bulkRNA, log2FC.bulkpro) %>%
  arrange(desc(avg_log2FC.allTumorcellsvsPT))
plotdata_df <- melt(plotdata_wide_df)
plotdata_df <- plotdata_df %>%
  mutate(data_type = ifelse(variable == "avg_log2FC.allTumorcellsvsPT", dataname_snrna,
                            ifelse(variable == "log2FC.bulkRNA", dataname_bulk_rna, dataname_bulk_protein))) %>%
  mutate(foldchange = 2^value) %>%
  mutate(y_plot = Gene)
summary(plotdata_df$foldchange)
plotdata_df <- plotdata_df %>%
  mutate(x_plot = ifelse(foldchange >= 10, 10, foldchange))
plotdata_df$y_plot <- factor(x = plotdata_df$Gene, levels = plotdata_wide_df$Gene)
plotdata_df$data_type <- factor(x = plotdata_df$data_type, levels = c(dataname_snrna, dataname_bulk_rna, dataname_bulk_protein))

## make colors
display.brewer.all()
colors_datatype <- brewer.pal(n = 4, name = "Set1")[c(1, 3, 4)]
names(colors_datatype) <- c(dataname_snrna, dataname_bulk_rna, dataname_bulk_protein)

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
p <- p + guides(fill = guide_legend(override.aes = list(size=4), nrow = 3, title = NULL, label.theme = element_text(size = 12)))
file2write <- paste0(dir_out, "Foldchanges", ".png")
png(file2write, width = 600, height = 900, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "Foldchanges", ".pdf")
pdf(file2write, width = 4.25, height = 7, useDingbats = F)
print(p)
dev.off()


