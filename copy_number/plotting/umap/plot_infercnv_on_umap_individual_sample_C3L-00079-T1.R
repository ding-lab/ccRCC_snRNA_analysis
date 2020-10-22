# Yige Wu @WashU Oct 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample (all cells)

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

# input seurat object master list -----------------------------------------
seurat_summary <- fread(input = "./Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  # filter(!(Aliquot %in% c("CPT0001540013", "CPT0002270013", "CPT0015810004", "CPT0023690004", "CPT0078510004", "CPT0086820004"))) %>%
  filter(!(Aliquot %in% c("CPT0075170013"))) %>%
  mutate(Path_seurat_object = paste0("./Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object


# input dependencies ------------------------------------------------------
## input CNV genes
ccrcc_cna_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input metadata
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## infercnv parameters
infercnv_run_id <- "20200305.v1"
dir_infercnv_all_runs <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
dir_infercnv_output <- paste0(dir_infercnv_all_runs, "Individual.", infercnv_run_id, "/")

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
# snRNA_aliquot_id_tmp <- "CPT0086820004"
for (snRNA_aliquot_id_tmp in "CPT0001260013") {
  ## get the readable aliquot id
  id_aliquot_wu <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## create output directory by aliquot
  dir_out1 <- paste0(dir_out, id_aliquot_wu, "/")
  dir.create(dir_out1)
  
  ## input individually processed seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  ## get umap coordates
  umap_tmp_df <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tmp_df$barcode <- rownames(umap_tmp_df)
  
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
  
  for (cytoband_tmp in "3p26.1") {
  # for (cytoband_tmp in unique(ccrcc_cna_genes_df$Cytoband)) {
    ## create output directory by chromosome region
    dir_out2 <- paste0(dir_out1, cytoband_tmp, "/")
    dir.create(dir_out2)
    
    ## run all the genes belonging to this chromosome region
    genes2plot <- ccrcc_cna_genes_df$Gene_Symbol[ccrcc_cna_genes_df$Cytoband == cytoband_tmp]
    
    for (gene_tmp in genes2plot) {
      if ((gene_tmp %in% cnv_state_df$V1)) {
        ## extract current gene related results from the infercnv result data frame
        infercnv_observe_gene_tab <- cnv_state_df %>%
          rename(gene_symbol = V1) %>%
          filter(gene_symbol == gene_tmp) %>%
          mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
          rename(copy_state = value) %>%
          select(gene_symbol, barcode, copy_state)
        
        ## add CNV state to the barcode - UMAP coordidate data frame
        tab2p <- umap_tmp_df
        tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
        
        ## map CNV state value to text
        tab2p$cnv_cat <- map_infercnv_state2category(copy_state = tab2p$copy_state)
        tab2p$cnv_cat %>% table()
        
        ## make cells with CNV appear on top
        tab2p <- tab2p %>%
          arrange(desc(cnv_cat))
        
        p <- ggplot() +
          geom_point(data = tab2p, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.5) +
          scale_color_manual(values = copy_number_colors)
        p <- p + theme_bw()
        p <- p + theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
        p <- p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                axis.text.y=element_blank(),axis.ticks=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank())
        p <- p + labs(color = paste0(gene_tmp, " Copy Number Status"))
        p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
        p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))
        file2write <- paste(dir_out2, id_aliquot_wu,".", gene_tmp, ".png", sep="")
        png(file2write, width = 800, height = 900, res = 150)
        print(p)
        dev.off()
        
        file2write <- paste0(dir_out2, id_aliquot_wu,".", gene_tmp, ".pdf")
        pdf(file2write, width = 8, height = 9, useDingbats = F)
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

