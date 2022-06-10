# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "ggrastr",
  "ggplot2",
  "readxl",
  "RColorBrewer",
  "ggrepel"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20210805.v1/MetaData_TumorCellOnlyReclustered.20210805.v1.tsv")
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## input known cnv genes
knowncnvgenes_df <- read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# pre-process --------------------------------------------------------------
## set genes to plot
genes2plot <- knowncnvgenes_df$Gene_Symbol
knowncnvgenes_df <- knowncnvgenes_df %>%
  mutate(chr_arm = paste0(Chromosome, Arm))
## set infercnv result direcotory
dir_infercnv_output <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
## set colors
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")
copy_number_colors <-  c("0 Copies" = PuBu_colors[9],
                         "1 Copy" = PuBu_colors[5],
                         "2 Copies" = "grey50",
                         "3 Copies" = PuRd_colors[5], 
                         "4 Copies" = PuRd_colors[7],
                         ">4 Copies" = PuRd_colors[9],
                         "Not Available" = "black")

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
for (infercnv_run_id in c("Individual.20200305.v1", "run.20210805")) {
# for (infercnv_run_id in c("Individual.20200305.v1")) {
  
  dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
  ## get aliquots to process
  aliquots2process <- list.files(dir_infercnv_run)
  aliquots2process <- aliquots2process[grepl(pattern = "CPT", x = aliquots2process)]
  
  # for (snRNA_aliquot_id_tmp in c("CPT0075130004")) {
    for (snRNA_aliquot_id_tmp in c("CPT0075120002", "CPT0075140002", "CPT0075130004")) {
    # for (snRNA_aliquot_id_tmp in aliquots2process) {
    ## get the case id for this aliquot to show in the title
    easy_id_tmp <- unique(barcode2umap_df$easy_id[barcode2umap_df$orig.ident == snRNA_aliquot_id_tmp])
    
    ## create output directory by aliquot
    dir_out1 <- paste0(dir_out, easy_id_tmp, "/")
    dir.create(dir_out1)
    
    scrublets_df <- barcode2scrublet_df %>%
      filter(Aliquot_WU == easy_id_tmp) %>%
      filter(predicted_doublet)
    barcodes_doublet <- scrublets_df$Barcode; length(barcodes_doublet)
    
    ## get umap coordates
    umap_df <- barcode2umap_df %>%
      filter(orig.ident == snRNA_aliquot_id_tmp) %>%
      filter(!(barcode_tumorcellreclustered %in% barcodes_doublet)) %>%
      mutate(barcode = barcode_tumorcellreclustered) %>%
      mutate(Text_TumorCluster = paste0("C", id_manual_cluster_w0+1))
    
    ## input infercnv CNV state results
    obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_run, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
    ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_run, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
    dim(obs_cnv_state_mat)
    dim(ref_cnv_state_mat)
    
    ## transform infercnv result wide data frame to long data frame
    cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
    rm(obs_cnv_state_mat)
    rm(ref_cnv_state_mat)
    
    for (gene_tmp in c("VHL", "SQSTM1")) {
    # for (gene_tmp in genes2plot) {
      chr_arm_tmp <- knowncnvgenes_df$chr_arm[knowncnvgenes_df$Gene_Symbol == gene_tmp]
      ## create output directory by chromosome region
      dir_out2 <- paste0(dir_out1, chr_arm_tmp, "/")
      dir.create(dir_out2, showWarnings = F)
      
      # file2write <- paste(dir_out2, easy_id_tmp, ".", gene_tmp, ".png", sep="")
      file2write <- paste(dir_out2, easy_id_tmp, ".", gene_tmp, ".pdf", sep="")
      if ((gene_tmp %in% cnv_state_df$V1) & !file.exists(file2write)) {
        ## extract current gene related results from the infercnv result data frame
        infercnv_observe_gene_tab <- cnv_state_df %>%
          rename(gene_symbol = V1) %>%
          filter(gene_symbol == gene_tmp) %>%
          mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
          rename(copy_state = value) %>%
          select(gene_symbol, barcode, copy_state)
        
        ## add CNV state to the barcode - UMAP coordidate data frame
        point_data_df <- merge(umap_df, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
        
        ## map CNV state value to text
        point_data_df$cnv_cat <- map_infercnv_state2category(copy_state = point_data_df$copy_state)
        point_data_df$cnv_cat %>% table()
        
        ## make cells with CNV appear on top
        point_data_df <- point_data_df %>%
          arrange(factor(cnv_cat, levels = c("2 Copies", ">4 Copies", "0 Copy", "1 Copy", "3 Copies")))
          # arrange(desc(cnv_cat))
        
        ## make text data
        cellnumber_percluster_df <- point_data_df %>%
          select(Text_TumorCluster) %>%
          table() %>%
          as.data.frame() %>%
          rename(Text_TumorCluster = ".")
        text_data_df <- point_data_df %>%
          filter(Text_TumorCluster %in% cellnumber_percluster_df$Text_TumorCluster[cellnumber_percluster_df$Freq >= 50]) %>%
          group_by(Text_TumorCluster) %>%
          summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        
        p <- ggplot() +
          geom_point_rast(data = point_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 0.8, size = 0.2) +
          scale_color_manual(values = copy_number_colors)
        p <- p + geom_text_repel(data = text_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, label = Text_TumorCluster))
        p <- p + theme_bw()
        p <- p + theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
        p <- p + theme(legend.position = "none")
        p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                axis.text.y=element_blank(),axis.ticks=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank())
        pdf(file2write, width = 2, height = 2, useDingbats = F)
        print(p)
        dev.off()
        
        # png(file2write, width = 800, height = 900, res = 150)
        # print(p)
        # dev.off()
      }
    }
  }
}


