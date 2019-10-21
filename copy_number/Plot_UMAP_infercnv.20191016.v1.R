# Yige Wu @WashU Oct 2019

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# functions ---------------------------------------------------------------
map_infercnv_state2category <- function(copy_state) {
  cnv_cat <- vector(mode = "character", length = length(copy_state))
  cnv_cat[is.na(copy_state)] <- "Not Available"
  cnv_cat[copy_state == 1] <- "Neutral"
  cnv_cat[copy_state == 0] <- "Complete Loss"
  cnv_cat[copy_state == 0.5] <- "Loss of one copy"
  cnv_cat[copy_state == 1.5] <- "Addition of one copy"
  cnv_cat[copy_state == 2] <- "Addition of two copies"
  return(cnv_cat)
}
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuBu_colors
RuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")

colors_manual <- c("Complete Loss" = PuBu_colors[9],
                   "Loss of one copy" = PuBu_colors[5], 
                   "Neutral" = PuBu_colors[3],
                   "Addition of one copy" = RuRd_colors[5],
                   "Addition of two copies" = RuRd_colors[9],
                   "Not Available" = "grey50")


# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
sample_id <- c("CPT0086350004")

# genes2plot <- c("VHL", "PBRM1", "BAP1", "SETD2", "HIF1A", "CA9", "JAK2")
genes2plot <- c("VHL", "PBRM1", "BAP1", "SETD2",
                "HIF1A",
                "CDKN2A", "PTEN", 
                "NEGR1",
                "QKI",
                "CADM2", 
                "PTPRD", "NRXN3",
                "PRKCI", 
                "MECOM",
                "MDM4",
                "MYC",
                "JAK2")

# input seurat object -----------------------------------------------------
object2plot <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")

# extract the UMAP coordinates for each cell ------------------------------
umap_tab <- FetchData(object2plot, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_tab$barcode_int <- rownames(umap_tab)
umap_tab <- umap_tab %>%
  filter(orig.ident %in% sample_id) %>%
  mutate(barcode = str_split_fixed(string = barcode_int, pattern = "_", n = 2)[,1])

# input infercnv observations ---------------------------------------------
tumor_cnv_state_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/inferCNV/outputs/integration.20190927.v1/CPT0086350004/infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt", data.table = F)
ref_cnv_state_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/inferCNV/outputs/integration.20190927.v1/CPT0086350004/infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.references.txt", data.table = F)
dim(tumor_cnv_state_mat)
dim(ref_cnv_state_mat)

cnv_state_df <- rbind(melt(tumor_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
table(cnv_state_df$value)

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
      scale_color_manual(values = colors_manual) +
      ggtitle(gene_tmp) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank())

    file2write <- paste(dir_out, sample_id,".", gene_tmp, ".FeaturePlot_CNA.", run_id, ".png", sep="")
    png(file2write, width = 1000, height = 800, res = 150)
    print(p)
    dev.off()
    
  } else {
    print(paste0(gene_tmp, " not in the infercnv result!"))
  }
  
}
