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
enricher_out <- readRDS(file = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/ora_msigdb_H_CP_tumorcells_up_degs/20210824.v1/ORA_Results.RDS")

# select non-overlapping pathways -----------------------------------------
View(enricher_out@termsim)
View(rownames(enricher_out@termsim))
## similarity value < 0.1, except for HALLMARK_GLYCOLYSIS
# pathways_selected <- c("HALLMARK_HYPOXIA", "HALLMARK_PROTEIN_SECRETION", "HALLMARK_GLYCOLYSIS", 
#                        "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", "WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR",
#                        "KEGG_RENAL_CELL_CARCINOMA", "REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_UV_RESPONSE_DN",
#                        "REACTOME_GLYCOGEN_METABOLISM", "KEGG_CELL_CYCLE", "HALLMARK_ADIPOGENESIS")
pathways_selected <- c("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "BIOCARTA_VITCB_PATHWAY", "KEGG_RNA_DEGRADATION", "REACTOME_SLC_TRANSPORTER_DISORDERS",
                       "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX",
                       "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", "WP_GLUCOCORTICOID_RECEPTOR_PATHWAY")

View(enricher_out@termsim[pathways_selected,])
View(enricher_out@termsim[pathways_selected,pathways_selected])

# make plot data ----------------------------------------------------------
plotdata_df <- enricher_out@result
plotdata_df <- plotdata_df %>%
  filter(Description %in% pathways_selected) %>%
  mutate(size_plot = Count) %>%
  mutate(x_plot = (size_plot/768)*100) %>%
  mutate(log10FDR = -log10(p.adjust))
# pathway_label_df <- data.frame(Description = pathways_selected,
#                                pathway_label = c("Hypoxia", "Protein secretion", "Glycolysis",
#                                                  "Antigen presentation of MHC I", "Primary cilium development",
#                                                  "Renal cell carcinoma signaling", "Pyruvate metabolism & TCA cycle", "Inteferon gamma response", "UV response (down-regulated)",
#                                                  "Glycogen metabolism", "Cell cycle", "Adipogenesis"))
pathway_label_df <- data.frame(Description = pathways_selected,
                               pathway_label = c("Hypoxia", "Protein secretion", "Glycolysis",
                                                 "Antigen presentation of MHC I", "Primary cilium development",
                                                 "Renal cell carcinoma signaling", "Pyruvate metabolism & TCA cycle", "Inteferon gamma response", "UV response (down-regulated)",
                                                 "Glycogen metabolism", "Cell cycle", "Adipogenesis"))

plotdata_df$y_plot <- mapvalues(x = plotdata_df$Description, from = pathway_label_df$Description, to = as.vector(pathway_label_df$pathway_label))
plotdata_df <- plotdata_df %>%
  arrange(x_plot)
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = plotdata_df$y_plot)

# plot enrichment map -----------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = log10FDR))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + guides(color = guide_legend(title = "-log10(FDR)", nrow = 2, title.position = "top"),
                size = guide_legend(title = "Gene count", nrow = 2, title.position = "top"))
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "bottom")
p
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 3, height = 3.5, useDingbats = F)
print(p)
dev.off()

# plot enrichment map -----------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + guides(size = guide_legend(title = "Gene count", nrow = 1, title.position = "top"))
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "bottom")
p
file2write <- paste(dir_out, "dotplot.colorlegend.pdf")
pdf(file2write, width = 3, height = 3.5, useDingbats = F)
print(p)
dev.off()


# plot enrichment map -----------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + guides(size = guide_legend(title = "Gene count", nrow = 1, title.position = "top"))
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "bottom")
p
file2write <- paste(dir_out, "dotplot.sizelegend.pdf")
pdf(file2write, width = 3, height = 3.5, useDingbats = F)
print(p)
dev.off()


