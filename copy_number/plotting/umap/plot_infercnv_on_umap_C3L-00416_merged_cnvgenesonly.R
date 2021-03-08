# Yige Wu @WashU March 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample (all cells)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat object master list -----------------------------------------
## input cnv data
cnv_state_united_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/run_infercnv/unite_infercnv_outputs14_cnvgenesonly/20201225.v1/infercnv.step14.outputs.20201225.v1.tsv")
## input CNV genes
ccrcc_cna_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input metadata
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# input seurat object, and get umap info -----------------------------------------------------
case_id <- "C3L-00416"
aliquot_show <- paste0(case_id, "_Merged")
srat <- readRDS(file = "./Data_Freezes/V1/snRNA/Merged_Seurat_Objects/C3L-00416.Tumor_Segments.Merged.20200319.v1.RDS")
umap_df <- FetchData(object = srat, vars = c("orig.ident", "UMAP_1", "UMAP_2"))
umap_df$barcode_integrated <- rownames(umap_df)
umap_df <- umap_df %>%
  mutate(barcode = str_split_fixed(string = barcode_integrated, pattern = "_", n = 2)[,1])
## filter cnv state
cnv_state_df <- cnv_state_united_df %>%
  filter(grepl(x = easy_id, pattern = case_id))

# make plot data ----------------------------------------------------------
for (cytoband_tmp in unique(ccrcc_cna_genes_df$Cytoband)) {
  ## create output directory by chromosome region
  dir_out1 <- paste0(dir_out, cytoband_tmp, "/")
  dir.create(dir_out1)
  
  ## run all the genes belonging to this chromosome region
  genes2plot <- ccrcc_cna_genes_df$Gene_Symbol[ccrcc_cna_genes_df$Cytoband == cytoband_tmp]
  
  for (gene_tmp in genes2plot) {
    if ((gene_tmp %in% cnv_state_df$gene_symbol)) {
      cnv_state_gene_df <- cnv_state_df %>%
        filter(gene_symbol == gene_tmp)
      
      ## add CNV state to the barcode - UMAP coordidate data frame
      plot_data_df <- umap_df
      plot_data_df <- merge(plot_data_df, cnv_state_gene_df, by.x = c("orig.ident", "barcode"), by.y = c("aliquot", "barcode"), all.x = T)
      plot_data_df$copy_state[is.na(plot_data_df$copy_state)] <- 1
      
      ## map CNV state value to text
      plot_data_df$cnv_cat <- map_infercnv_state2category(copy_state = plot_data_df$copy_state)
      plot_data_df$cnv_cat %>% table()
      
      ## make cells with CNV appear on top
      plot_data_df <- rbind(plot_data_df[plot_data_df$copy_state == 1,], plot_data_df[plot_data_df$copy_state != 1,])
      
      p <- ggplot() +
        geom_point_rast(data = plot_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.3) +
        scale_color_manual(values = copy_number_colors)
      p <- p + theme_bw()
      p <- p + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
      p <- p + theme(legend.position = "none")
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      file2write <- paste(dir_out1, aliquot_show,".", gene_tmp, ".nolegend", ".pdf", sep="")
      pdf(file2write, width = 4, height = 4, useDingbats = F)
      print(p)
      dev.off()
    }
    if (!(gene_tmp %in% cnv_state_df$gene_symbol)) {
      sink(file = paste0(dir_out1, "genes_not_in_infercnv_output.txt"))
      cat(paste0(gene_tmp, " not in the infercnv result!\n"))
      sink()
    }
  }
}

