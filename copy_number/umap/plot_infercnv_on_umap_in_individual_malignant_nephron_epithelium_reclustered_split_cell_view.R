# Yige Wu @WashU Feb 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
### filter out the aliquots without CNV result yet
srat_paths <- srat_paths %>%
  filter(!(Aliquot %in% c("CPT0001540013", "CPT0002270013", "CPT0015810004", "CPT0078510004", "CPT0086820004", "CPT0079410004")))
## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
infercnv_run_id <- "20200207.v1"
dir_infercnv_all_runs <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
dir_infercnv_output <- paste0(dir_infercnv_all_runs, "Individual.", infercnv_run_id, "/")

for (snRNA_aliquot_id_tmp in srat_paths$Aliquot) {
  dir_out1 <- paste0(dir_out, snRNA_aliquot_id_tmp, "/")
  
  if (dir.exists(dir_out1)) {
    next()
  }
  ## get case id
  case_id_tmp <- meta_tab$Case[meta_tab$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## create output directory by aliquot
  dir.create(dir_out1)
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  ## get umap coordates
  umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  ## get labels associated with umap coordinates
  p <- DimPlot(seurat_object, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  ## input infercnv CNV state results
  obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(obs_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  ## transform infercnv result wide data frame to long data frame
  cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  rm(obs_cnv_state_mat)
  rm(ref_cnv_state_mat)
  
  for (chr_region_tmp in unique(ccrcc_cna_genes_df$chr_region)) {
    ## create output directory by chromosome region
    dir_out2 <- paste0(dir_out1, chr_region_tmp, "/")
    dir.create(dir_out2)
    
    ## run all the genes belonging to this chromosome region
    genes2plot <- ccrcc_cna_genes_df$gene_symbol[ccrcc_cna_genes_df$chr_region == chr_region_tmp]
    
    for (gene_tmp in genes2plot) {
      file2write <- paste(dir_out2, snRNA_aliquot_id_tmp, ".nephron_epithelium_reclustered.", "20200217.v1",
                          ".umap.", gene_tmp, ".", run_id, ".png", sep="")
      if ((gene_tmp %in% cnv_state_df$V1) & !file.exists(file2write)) {
        ## extract current gene related results from the infercnv result data frame
        infercnv_observe_gene_tab <- cnv_state_df %>%
          rename(gene_symbol = V1) %>%
          filter(gene_symbol == gene_tmp) %>%
          mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
          rename(copy_state = value) %>%
          select(gene_symbol, barcode, copy_state)
        
        ## add gene level CNV state to the barcode - UMAP coordidate data frame
        tab2p <- umap_tab
        tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
        
        ## map CNV state value to text
        tab2p$cnv_cat <- map_infercnv_state2category(copy_state = tab2p$copy_state)
        tab2p$cnv_cat %>% table()
        
        ## make cells with CNV appear on top
        tab2p <- tab2p %>%
          arrange(desc(cnv_cat))
        
        ## get the case id for this aliquot to show in the title
        p <- ggplot() +
          geom_point(data = tab2p, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.3) +
          scale_color_manual(values = copy_number_colors)
        p <- p + facet_grid(.~cnv_cat, drop = T)
        p <- p + ggtitle(paste0("Case: ", case_id_tmp, "   Aliquot: ",  snRNA_aliquot_id_tmp), 
                         subtitle = paste0(gene_tmp, " Copy Number Status"))
        p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
        p <- p + theme_bw()
        p <- p + theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
        p <- p + theme(legend.position = "top")
        p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                axis.text.y=element_blank(),axis.ticks=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank())
        if (length(unique(tab2p$cnv_cat)) == 2) {
          png(file2write, width = 1500, height = 900, res = 150)
        } else if (length(unique(tab2p$cnv_cat)) == 3) {
          png(file2write, width = 2000, height = 900, res = 150)
        } else {
          png(file2write, width = 2800, height = 900, res = 150)
        }
        print(p)
        dev.off()
      }
      if (!(gene_tmp %in% cnv_state_df$V1)) {
        sink(file = paste0(dir_out2, "genes_not_in_infercnv_output.txt"))
        cat(paste0(gene_tmp, " not in the infercnv result!\n"))
        sink()
      }
    }
  }
}


