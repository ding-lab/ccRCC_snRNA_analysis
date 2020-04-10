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
integration_id <- "20191021.v1"
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"

# set aliquot id ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# set genes to plot -------------------------------------------------------
genes2plot <- ccrcc_cna_genes_df$gene_symbol
genes2plot <- "VHL"
genes2plot <- "HIF1A"
genes2plot <- c("VHL", "HIF1A", "EPAS1")
genes2plot <- c("PBRM1", "BAP1", "SETD2")

# plot by each aliquot ----------------------------------------------------
# snRNA_aliquot_id_tmp <- c("CPT0019130004")
# snRNA_aliquot_id_tmp <- c("CPT0025890002")
# snRNA_aliquot_id_tmp <- c("CPT0086350004")
# snRNA_aliquot_id_tmp <- c("CPT0010110013")
# snRNA_aliquot_id_tmp <- c("CPT0001180011")
# snRNA_aliquot_id_tmp <- c("CPT0001260013")
# snRNA_aliquot_id_tmp <- c("CPT0075140002")
snRNA_aliquot_id_tmp <- c("CPT0001220012")
snRNA_aliquot_id_tmp <- "CPT0001220012"

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # input individually processed seurat object ---------------------------------------------
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  # get umap coordinates ----------------------------------------------------
  umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  # input infercnv observations ---------------------------------------------
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
      
      tab2p <- tab2p %>%
        arrange(desc(cnv_cat))
      
      # plot a UMAP plot for copy number metric ---------------------------------
      p <- ggplot(tab2p, aes(UMAP_1, UMAP_2)) +
        geom_point(aes(color=cnv_cat), alpha = 1, size = 0.3) +
        scale_color_manual(values = copy_number_colors) +
        ggtitle(paste0(snRNA_aliquot_id_tmp, "_", gene_tmp, "_Copy_Number_Status"))
      p <- p + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      file2write <- paste(dir_out, snRNA_aliquot_id_tmp,"_Individual_Clustered_", gene_tmp, ".FeaturePlot_CNA.", run_id, ".png", sep="")
      png(file2write, width = 1100, height = 800, res = 150)
      print(p)
      dev.off()
    } else {
      print(paste0(gene_tmp, " not in the infercnv result!"))
    }
  }
}

