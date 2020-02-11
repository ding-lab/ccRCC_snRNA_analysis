# Yige Wu @WashU Oct 2019

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat object -----------------------------------------------------
object2plot <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")

# set aliquot id ----------------------------------------------------------
# snRNA_aliquot_id_tmp <- c("CPT0019130004")
# snRNA_aliquot_id_tmp <- c("CPT0025890002")
# snRNA_aliquot_id_tmp <- c("CPT0086350004")
# snRNA_aliquot_id_tmp <- c("CPT0010110013")
# snRNA_aliquot_id_tmp <- c("CPT0001180011")
snRNA_aliquot_id_tmp <- c("CPT0001260013")
# snRNA_aliquot_id_tmp <- c("CPT0075140002")

# # input infercnv predicted CNA state by subcluster by gene ---------------------------------------------
# infercnv_subcluster_cna_state_by_gene <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.20191021.v1/", snRNA_aliquot_id_tmp, "/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"), data.table = F)
# infercnv_subcluster_cna_state_by_region <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.20191021.v1/", snRNA_aliquot_id_tmp, "/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"), data.table = F)
# infercnv_cna_state_by_region <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.20191021.v1/", snRNA_aliquot_id_tmp, "/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"), data.table = F)
# 
# # search for subcluster specific CNA --------------------------------------
# infercnv_subcluster_cna_state_filtered <- infercnv_subcluster_cna_state %>%
#   filter(gene %in% ccrcc_cna_genes_df$gene_symbol)
# 
# infercnv_subcluster_cna_state_filtered_mat <- dcast(data = infercnv_subcluster_cna_state_filtered, formula = gene ~ cell_group_name, value.var = "state")

# extract the UMAP coordinates for each cell ------------------------------
umap_tab <- FetchData(object2plot, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_tab$barcode_int <- rownames(umap_tab)
umap_tab <- umap_tab %>%
  filter(orig.ident %in% snRNA_aliquot_id_tmp) %>%
  mutate(barcode = str_split_fixed(string = barcode_int, pattern = "_", n = 2)[,1])

# plot subcluster groupings on UMAP ---------------------------------------
tab2p <- umap_tab
tab2p <- merge(tab2p, infercnv_subcluster_groupings, by.x = c("barcode"), by.y = c("cell"), all.x = T)

# # plot a UMAP plot for copy number metric ---------------------------------
# p <- ggplot(tab2p, aes(UMAP_1, UMAP_2)) +
#   geom_point(aes(color=cell_group_name), alpha = 0.7) +
#   # scale_color_manual(values = colors_manual) +
#   ggtitle(paste0(snRNA_aliquot_id_tmp, "_", "Subclusters")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_blank()) +
#   guides(color = FALSE)
# 
# file2write <- paste(dir_out, snRNA_aliquot_id_tmp, "DimPlot_Subcluster_Groupings.", run_id, ".png", sep="")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()


# set genes to plot -------------------------------------------------------
genes2plot <- c("HIF1A")
genes2plot <- c("VHL", "PBRM1", "BAP1", "SETD2", "SQSTM1", "HIF1A")


# input infercnv observations ---------------------------------------------
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
integration_id <- "20191021.v1"
tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
dim(tumor_cnv_state_mat)
dim(ref_cnv_state_mat)

cnv_state_df <- rbind(melt(tumor_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))

for (gene_tmp in genes2plot) {
  if (gene_tmp %in% cnv_state_df$V1) {
    infercnv_observe_gene_tab <- cnv_state_df %>%
      rename(gene_symbol = V1) %>%
      filter(gene_symbol == gene_tmp) %>%
      mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
      rename(copy_state = value) %>%
      select(gene_symbol, barcode, copy_state)
    
    tab2p <- umap_tab
    tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
    tab2p$cnv_cat <- map_infercnv_state2category(copy_state = tab2p$copy_state)
    tab2p$cnv_cat %>% table()
    
    # plot a UMAP plot for copy number metric ---------------------------------
    p <- ggplot(tab2p, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color=cnv_cat), alpha = 0.7) +
      scale_color_manual(values = copy_number_colors) +
      ggtitle(paste0(snRNA_aliquot_id_tmp, "_", gene_tmp)) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank())

    file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".", gene_tmp, ".FeaturePlot_CNA.", run_id, ".png", sep="")
    png(file2write, width = 1000, height = 800, res = 150)
    print(p)
    dev.off()
    
  } else {
    print(paste0(gene_tmp, " not in the infercnv result!"))
  }
}
