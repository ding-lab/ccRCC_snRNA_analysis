# Yige Wu @WashU Jun 2021

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
## input the enricher object
# enricher_out <- readRDS(file = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/ora_msigdb_H_CP_tumorcells_up_degs/20210629.v1/ORA_Results.RDS")
enricher_out <- readRDS(file = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/ora_msigdb_H_CP_tumorcells_down_degs_w_background/20210825.v1/ORA_Results.RDS")
enricher_out_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/ora_msigdb_H_CP_tumorcells_down_degs_w_background/20210825.v1/ORA_Results.tsv")

# select non-overlapping pathways -----------------------------------------
View(enricher_out@termsim)
View(rownames(enricher_out@termsim))
## similarity value < 0.1
pathways_selected <- c("REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT", "WP_EICOSANOID_METABOLISM_VIA_LIPO_OXYGENASES_LOX", "NABA_ECM_AFFILIATED", "WP_ZINC_HOMEOSTASIS", "HALLMARK_ESTROGEN_RESPONSE_LATE")
pathways_examine <- intersect(enricher_out_df$ID, rownames(enricher_out@termsim))
View(t(enricher_out@termsim[pathways_selected,pathways_examine]))
View(enricher_out@termsim[pathways_selected,pathways_selected])

# make plot data ----------------------------------------------------------
plotdata_df <- enricher_out_df
plotdata_df <- plotdata_df %>%
  filter(Description %in% pathways_selected) %>%
  mutate(size_plot = Count) %>%
  mutate(x_plot = (size_plot/768)*100) %>%
  mutate(log10FDR = -log10(p.adjust))
pathway_label_df <- data.frame(Description = pathways_selected,
                               pathway_label = c("SLC transmembrane trasnport", "Eicosanoid metabolism", "Extracellular matrix proteins", "Zinc homeostasis", "Late response to estrogen"))

plotdata_df$y_plot <- mapvalues(x = plotdata_df$Description, from = pathway_label_df$Description, to = as.vector(pathway_label_df$pathway_label))
plotdata_df <- plotdata_df %>%
  arrange(x_plot)
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = plotdata_df$y_plot)

# plot enrichment map -----------------------------------------------------
## make colors for log10FDR
colors_fdr <- brewer.pal(n = 11, name = "Spectral")
colors_fdr <- rev(colors_fdr[10, 6, 5, 4, 3, 2, 1])
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = log10FDR))
# p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"), breaks = c(2, 4, 6), guide = guide_colourbar(direction = "horizontal", title = NULL))
p <- p + scale_size_continuous(limits = c(0, 40), breaks = c(10, 20, 30, 40), guide = guide_legend(direction = "horizontal", title = NULL, nrow = 2, byrow = T))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
# p <- p + guides(color = guide_legend(title = "-log10(FDR)", title.position = "top"),
#                 size = guide_legend(title = "Gene count", nrow = 1, title.position = "top"))
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "right", legend.box = "vertical")
p
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 4.5, height = 1.4, useDingbats = F)
print(p)
dev.off()

# # plot enrichment map -----------------------------------------------------
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot))
# p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
# p <- p + theme_light(base_size = 12)
# p <- p + xlab(label = "Gene ratio (%)")
# p <- p + guides(size = guide_legend(title = "Gene count", nrow = 1, title.position = "top"))
# p <- p + theme(axis.text = element_text(color = "black"),
#                axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
#                legend.position = "bottom")
# p
# file2write <- paste(dir_out, "dotplot.colorlegend.pdf")
# pdf(file2write, width = 3, height = 3.5, useDingbats = F)
# print(p)
# dev.off()
# 
# 
# # plot enrichment map -----------------------------------------------------
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot))
# p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
# p <- p + theme_light(base_size = 12)
# p <- p + xlab(label = "Gene ratio (%)")
# p <- p + guides(size = guide_legend(title = "Gene count", nrow = 1, title.position = "top"))
# p <- p + theme(axis.text = element_text(color = "black"),
#                axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
#                legend.position = "bottom")
# p
# file2write <- paste(dir_out, "dotplot.sizelegend.pdf")
# pdf(file2write, width = 3, height = 3.5, useDingbats = F)
# print(p)
# dev.off()
# 
# 
