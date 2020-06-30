# Yige Wu @WashU Feb 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(readxl)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for individual seurat object
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## set infercnv result direcotory
infercnv_run_id <- "20200305.v1"
dir_infercnv_all_runs <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
dir_infercnv_output <- paste0(dir_infercnv_all_runs, "Individual.", infercnv_run_id, "/")
## input barcode-to-manual grouping 
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv", data.table = F)
## input known cnv genes
knowncnvgenes_df <- read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## set genes to plot
genes2plot <- knowncnvgenes_df$Gene_Symbol

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
# for (snRNA_aliquot_id_tmp in "CPT0075130004") {
for (snRNA_aliquot_id_tmp in srat_paths$Aliquot) {
  ## get the case id for this aliquot to show in the title
  id_aliquot_wu_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## create output directory by aliquot
  dir_out1 <- paste0(dir_out, id_aliquot_wu_tmp, "/")
  dir.create(dir_out1)
  
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  ## get barcode to cnv info for this aliquot
  barcode2manualsubcluster_aliquot_df <- barcode2manualsubcluster_df %>%
    filter(orig.ident == snRNA_aliquot_id_tmp) %>%
    mutate(Name_Cluster = paste0("C", (Id_TumorManualCluster+1)))
  unique(barcode2manualsubcluster_aliquot_df$Name_Cluster)
  
  ## get umap coordates
  umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  umap_tab$Name_Cluster <- mapvalues(x = umap_tab$barcode, from = barcode2manualsubcluster_aliquot_df$individual_barcode, to = as.vector(barcode2manualsubcluster_aliquot_df$Name_Cluster))
  
  ## get labels associated with umap coordinates
  p <- DimPlot(seurat_object, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  label_data$Name_Cluster <- mapvalues(x = label_data$ident, from = umap_tab$ident, to = as.vector(umap_tab$Name_Cluster))
  label_data <- label_data[!duplicated(label_data$Name_Cluster),]
  
  ## input infercnv CNV state results
  obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(obs_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  ## transform infercnv result wide data frame to long data frame
  cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  rm(obs_cnv_state_mat)
  rm(ref_cnv_state_mat)
  
  for (gene_tmp in genes2plot) {
    file2write <- paste(dir_out1, snRNA_aliquot_id_tmp, ".nephron_epithelium_reclustered.", "20200217.v1",
                        ".umap.", gene_tmp, ".", run_id, ".png", sep="")
    if ((gene_tmp %in% cnv_state_df$V1) & !file.exists(file2write)) {
      ## extract current gene related results from the infercnv result data frame
      infercnv_observe_gene_tab <- cnv_state_df %>%
        rename(gene_symbol = V1) %>%
        filter(gene_symbol == gene_tmp) %>%
        mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
        rename(copy_state = value) %>%
        select(gene_symbol, barcode, copy_state)
      
      ## add CNV state to the barcode - UMAP coordidate data frame
      tab2p <- umap_tab
      tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
      
      ## map CNV state value to text
      tab2p$cnv_cat <- map_infercnv_state2category(copy_state = tab2p$copy_state)
      tab2p$cnv_cat %>% table()
      
      ## make cells with CNV appear on top
      tab2p <- tab2p %>%
        arrange(desc(cnv_cat))
      
      
      p <- ggplot() +
        geom_point(data = tab2p, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.3) +
        scale_color_manual(values = copy_number_colors)
      p <- p + ggtitle(paste0(id_aliquot_wu_tmp, " Tumor-Cell-Only Clusters"), 
                       subtitle = paste0(gene_tmp, " Copy Number Status"))
      p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = Name_Cluster))
      p <- p + theme_bw()
      p <- p + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
      p <- p + theme(legend.position = "top")
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      png(file2write, width = 800, height = 900, res = 150)
      print(p)
      dev.off()
    }
  }
}


