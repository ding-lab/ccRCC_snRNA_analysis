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
## input gene-to-pathway
# gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pathway/unite_PBRM1_BAP1_vs_NonMutants_DAP_ORA/20210628.v2/PBRM1_BAP1_vs_NonMutants.DAPGene2TopPathway.20210628.v2.tsv")
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pathway/unite_PBRM1_BAP1_vs_NonMutants_DAP_ORA_28samples/20211011.v1/PBRM1_BAP1_vs_NonMutants.DAPGene2TopPathway.20211011.v1.tsv")

## input pbrm1 daps
daps_pbrm1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/annotate_peaks/annotate_pbrm1_vs_nonmutant_daps_28samples/20211011.v1/PBRM1_DAP2Gene.EnhancerPromoter.20211011.v1.tsv")
daps_bap1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_peaks/annotate_BAP1_vs_NonMutant_daps_28samples/20211011.v1/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20211011.v1.tsv")
## input pathway enrichment results
ora_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pathway/unite_PBRM1_BAP1_vs_NonMutants_DAP_ORA_28samples/20211011.v1/PBRM1_BAP1_vs_NonMutants.DAPTopPathways2Genes.20211011.v1.tsv")

# make plot data ----------------------------------------------------------
colnames_merge <- intersect(colnames(daps_pbrm1_df), colnames(daps_bap1_df))
plotdata_source_df <- rbind(daps_pbrm1_df[,colnames_merge] %>%
                              mutate(DAP_group = paste0("PBRM1_", DAP_direction)),
                            daps_bap1_df[,colnames_merge] %>%
                              mutate(DAP_group = paste0("BAP1_", DAP_direction)))
plotdata_source_df <- merge(x = plotdata_source_df, y = gene2pathway_df, by.x = c("Gene"), by.y = c("GeneSymbol"), all.x = T)  
## 
plotdata_df <- plotdata_source_df %>%
  filter(!is.na(GeneSet_Name)) %>%
  group_by(DAP_group) %>%
  select(GeneSet_Name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)
plotdata_df <- merge(x = plotdata_df, 
                     y = ora_df %>%
                       select(Description, p.adjust, Comparison, rank_pvalue),
                     by.x = c("GeneSet_Name", "DAP_group"), by.y = c("Description", "Comparison"), all.x = T)
plotdata_df <- plotdata_df %>%
  mutate(log10FDR = -log10(p.adjust))
## order the pathways
pathway_selected <- c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                      "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                      "WP_FOCAL_ADHESION", "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA")
pathway_ordered <- unique(pathway_selected)
plotdata_df$GeneSet_Name <- factor(x = plotdata_df$GeneSet_Name, levels = rev(pathway_ordered))
## make color patette
color_pal <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
## renmake pathway labels
pathway_label_df <- data.frame(label = c("mTORC1 signaling", "EGF-EGFR signaling", "Rho GTPase cycle", "DNA repair", 
                                         "Focal adhesion", "TGFA signaling via NFKB", "EPHA forward signaling", "TGFB signaling for EMT", 
                                         "Rac1 GTPase cycle", "Hypoxia"),
                               pathway_name = c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                                                "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA"))

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq, color = rank_pvalue), shape = 16)
p <- p + scale_size_area()
p <- p + scale_color_gradientn(colours = rev(color_pal), breaks = c(1, 3, 5, 7, 9))
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, color = "black"),
               axis.title = element_blank(),
               axis.text.y = element_text(color = "black"),
               legend.position = "right", legend.direction = "horizontal")
p <- p + guides(color = guide_legend(title = "Rank by\np-value", nrow = 2, title.position = "top"),
                size = guide_legend(title = "Number of\npeaks", nrow = 2, title.position = "top"))
p <- p + scale_y_discrete(breaks = pathway_label_df$pathway_name, labels = pathway_label_df$label)
file2write <- paste0(dir_out, "bubble.colorByRankpvalue.legend.pdf")
pdf(file2write, width = 3.85, height = 2.75, useDingbats = F)
print(p)
dev.off()


# plot by log10FDR ----------------------------------------------------------------
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq))
# p <- p + scale_size_area()
# file2write <- paste0(dir_out, "bubble.plain.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
# 
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq, fill = log10FDR), color = "black", shape = 21)
# p <- p + scale_size_area()
# # p <- p + scale_color_gradient(low = "white smoke", high = "red")
# p <- p + scale_fill_gradientn(colours = color_pal)
# p <- p + theme_classic()
# p
# file2write <- paste0(dir_out, "bubble.colorBylog10FDR.pdf")
# pdf(file2write, width = 6, height = 4, useDingbats = F)
# print(p)
# dev.off()
# file2write <- paste0(dir_out, "bubble.colorlog10FDR.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()

