# Yige Wu @WashU Oct 2019

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set infercnv input id and directory --------------------------------------------
integration_id <- "20191021.v1"
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"


# set genes to plot -------------------------------------------------------
cna_gene1 <- "VHL"
cna_gene2 <- "PBRM1"

# set aliquot id ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# set recluster_tumor_id --------------------------------------------------
recluster_tumor_id <- "20191119.v1"

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/recluster_tumor/", recluster_tumor_id, "/", Aliquot, FACS, 
                                     "/", Aliquot, FACS, ".Malignant_Reclustered.", recluster_tumor_id, ".RDS"))
seurat_summary2process$Path_seurat_object


# set genes to plot -------------------------------------------------------
genes2plot <- ccrcc_cna_genes_df$gene_symbol

# plot by each aliquot ----------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0086350004"

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  dir_out_tmp <- paste0(dir_out, snRNA_aliquot_id_tmp, "/")
  dir.create(dir_out_tmp)
  
  ## input individually clustered tumor cells
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  ## get umap cell coordinates
  umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  ## get umap label coordinates
  p <- DimPlot(seurat_object, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  # make the subplot for CNA mapping -----------------------------------
  
  ## input infercnv observations
  tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(tumor_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  cnv_state_df <- rbind(melt(tumor_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  cnv_state_df %>% head()
  
  ## filter by cna genes
  cnv_state_df <- cnv_state_df %>%
    filter(V1 %in% c(cna_gene1, cna_gene2)) %>%
    mutate(cna_group = ifelse(value == 1, "Neutral",
                              ifelse(value < 1, "Loss", "Gain")))
  
  gene1_cnv_state_df <- cnv_state_df %>%
    filter(V1 == cna_gene1)
  
  gene2_cnv_state_df <- cnv_state_df %>%
    filter(V1 == cna_gene2)
  
  ## get desired CNA type of gene1 and gene2
  desired_gene1_cna <- ccrcc_cna_genes_df$gene_cna_type[ccrcc_cna_genes_df$gene_symbol == cna_gene1]
  desired_gene2_cna <- ccrcc_cna_genes_df$gene_cna_type[ccrcc_cna_genes_df$gene_symbol == cna_gene2]
  
  ## categorize cells by CNA state in cna_gene1 and cna_gene2
  tab2p1 <- umap_tab
  tab2p1$gene1_cna_group <- mapvalues(x = tab2p1$barcode, from = gene1_cnv_state_df$variable, to = gene1_cnv_state_df$cna_group)
  tab2p1$gene2_cna_group <- mapvalues(x = tab2p1$barcode, from = gene2_cnv_state_df$variable, to = gene2_cnv_state_df$cna_group)
  
  tab2p1 <- tab2p1 %>%
    mutate(combined_cna_group = ifelse(gene1_cna_group == "Neutral" & gene2_cna_group == "Neutral", paste0(cna_gene1, "_Neutral", "&", cna_gene2, "_Neutral"),
                                       ifelse(gene1_cna_group == desired_gene1_cna & gene2_cna_group == desired_gene2_cna, paste0(cna_gene1, "_", desired_gene1_cna, "&", cna_gene2, "_", desired_gene2_cna),
                                              ifelse(gene1_cna_group == desired_gene1_cna & gene2_cna_group == "Neutral", paste0(cna_gene1, "_", desired_gene1_cna, " Only"),
                                                     ifelse(gene1_cna_group == "Neutral" & gene2_cna_group == desired_gene2_cna, paste0(cna_gene2, "_", desired_gene2_cna, " Only"), "Other")))))
  
  ## make CNA subplot
  combined_cna_colors <- c("grey70", "#7b3294", "#018571", "#a6611a", "#e78ac3")
  names(combined_cna_colors) <- c(paste0(cna_gene1, "_Neutral", "&", cna_gene2, "_Neutral"),
                                  paste0(cna_gene1, "_", desired_gene1_cna, "&", cna_gene2, "_", desired_gene2_cna),
                                  paste0(cna_gene1, "_", desired_gene1_cna, " Only"),
                                  paste0(cna_gene2, "_", desired_gene2_cna, " Only"),
                                  "Other")
  
  p1 <- ggplot() +
    geom_point(data = tab2p1, mapping = aes(UMAP_1, UMAP_2, color=combined_cna_group), alpha = 0.7, size = 1) +
    scale_color_manual(values = combined_cna_colors) +
    ggtitle(paste0(snRNA_aliquot_id_tmp, "_", cna_gene1, "&", cna_gene2, "_Copy_Number_Status"))
  p1 <- p1 + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
  p1 <- p1 + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p1 <- p1 + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
  file2write <- paste(dir_out_tmp, snRNA_aliquot_id_tmp,".Malignant_Reclustered.CNA_", cna_gene1, "_", cna_gene2,  ".combined.", run_id, ".png", sep="")
  png(file2write, width = 1200, height = 800, res = 150)
  print(p1)
  dev.off()
  
  p1_split <- ggplot() +
    geom_point(data = tab2p1, mapping = aes(UMAP_1, UMAP_2, color=combined_cna_group), alpha = 0.7, size = 1) +
    scale_color_manual(values = combined_cna_colors) +
    ggtitle(paste0(snRNA_aliquot_id_tmp, "_", cna_gene1, "&", cna_gene2, "_Copy_Number_Status"))
  p1_split <- p1_split + facet_grid(~combined_cna_group)
  p1_split <- p1_split + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
  p1_split <- p1_split + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p1_split <- p1_split + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                           axis.text.y=element_blank(),axis.ticks=element_blank(),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),
                           legend.position = "none")
  file2write <- paste(dir_out_tmp, snRNA_aliquot_id_tmp,".Malignant_Reclustered.CNA_", cna_gene1, "_", cna_gene2,  ".split.", run_id, ".png", sep="")
  png(file2write, width = 3500, height = 800, res = 150)
  print(p1_split)
  dev.off()

  # make the subplot for mutation mapping -----------------------------------
  mutation_map_tab <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/10Xmapping/outputs/", snRNA_aliquot_id_tmp, "/", snRNA_aliquot_id_tmp, "_mapping_heatmap_0.txt"), data.table = F)
  mutation_map_tab <- mutation_map_tab %>%
    mutate(gene_symbol = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,1])
  mutation_map_tab.m <- melt(mutation_map_tab, id.vars = c("gene_symbol", "Mutatation"))
  mutation_map_tab.m <- mutation_map_tab.m %>%
    mutate(allele_type = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,3])
  
  mutation_map_tab.var <- mutation_map_tab.m %>%
    filter(allele_type == "Var") %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable) %>%
    mutate(mutation = paste0(str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,1], "-", str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,2]))
  
  mut_genes <- mutation_map_tab.var %>%
    select(gene_symbol) %>%
    unique() 
  mut_genes <- mut_genes$gene_symbol
  
  mutation_map_tab.ref <- mutation_map_tab.m %>%
    filter(allele_type == "Ref") %>%
    filter(Mutatation %in% mutation_map_tab.var$Mutatation) %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable)
  
  tab2p2 <- umap_tab
  tab2p2$read_type <- mapvalues(x = tab2p2$barcode, from = mutation_map_tab.var$barcode, to = mutation_map_tab.var$mutation)
  tab2p2$read_type[tab2p2$read_type %in% tab2p2$barcode & tab2p2$barcode %in% mutation_map_tab.ref$barcode] <- "Ref"
  tab2p2$read_type[tab2p2$read_type %in% tab2p2$barcode] <- "NA"
  table(tab2p2$read_type)
  tab2p2$combined_cna_group <- mapvalues(x = tab2p2$barcode, from = tab2p1$barcode, to = tab2p1$combined_cna_group)
  
  mut_label_data <- tab2p2 %>%
    filter(!(read_type %in% c("NA", "Ref")))
  
  colors_read_type <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(mutation_map_tab.var$mutation))], "green", "grey70")
  names(colors_read_type) <- c(unique(mutation_map_tab.var$mutation), "Ref", "NA")
  
  p2 <- ggplot() +
    geom_point(data = tab2p2[tab2p2$read_type == "NA",], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 0.3) +
    geom_point(data = tab2p2[tab2p2$read_type == "Ref",], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 0.3) +
    geom_point(data = tab2p2[!(tab2p2$read_type %in% c("NA", "Ref")),], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 2) +
    scale_color_manual(values = colors_read_type) +
    ggtitle(paste0(snRNA_aliquot_id_tmp, "_Mutation_Status")) +
    theme_bw()
  p2 <- p2 + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
  p2 <- p2 + geom_text_repel(data = mut_label_data, mapping = aes(UMAP_1, UMAP_2, label = read_type), box.padding = 0.5,  force = 2)
  p2 <- p2 + facet_grid(~combined_cna_group)
  p2 <- p2 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line=element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.position = "none")
  file2write <- paste(dir_out_tmp, snRNA_aliquot_id_tmp,".Malignant_Reclustered.", "Mutation_All.", run_id, ".png", sep="")
  png(file2write, width = 3500, height = 800, res = 150)
  print(p2)
  dev.off()
  
  for (mut_gene in mut_genes) {
    dir_out_tmp2 <- paste0(dir_out_tmp, "Mutations/")
    dir.create(dir_out_tmp2)
    
    tab2p2 <- umap_tab
    tab2p2$read_type <- mapvalues(x = tab2p2$barcode, from = mutation_map_tab.var$barcode[mutation_map_tab.var$gene_symbol == mut_gene], to = mutation_map_tab.var$mutation[mutation_map_tab.var$gene_symbol == mut_gene])
    tab2p2$read_type[tab2p2$read_type %in% tab2p2$barcode & tab2p2$barcode %in% mutation_map_tab.ref$barcode[mutation_map_tab.var$gene_symbol == mut_gene]] <- "Ref"
    tab2p2$read_type[tab2p2$read_type %in% tab2p2$barcode] <- "NA"
    table(tab2p2$read_type)
    tab2p2$combined_cna_group <- mapvalues(x = tab2p2$barcode, from = tab2p1$barcode, to = tab2p1$combined_cna_group)
    
    mut_label_data <- tab2p2 %>%
      filter(!(read_type %in% c("NA", "Ref")))
    
    colors_read_type <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(mutation_map_tab.var$mutation[mutation_map_tab.var$gene_symbol == mut_gene]))], "green", "grey70")
    names(colors_read_type) <- c(unique(mutation_map_tab.var$mutation[mutation_map_tab.var$gene_symbol == mut_gene]), "Ref", "NA")
    
    p2 <- ggplot() +
      geom_point(data = tab2p2[tab2p2$read_type == "NA",], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 0.3) +
      geom_point(data = tab2p2[tab2p2$read_type == "Ref",], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 0.3) +
      geom_point(data = tab2p2[!(tab2p2$read_type %in% c("NA", "Ref")),], aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 1, size = 2) +
      scale_color_manual(values = colors_read_type) +
      ggtitle(paste0(snRNA_aliquot_id_tmp, "_Mutation_Status")) +
      theme_bw()
    p2 <- p2 + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
    p2 <- p2 + geom_text_repel(data = mut_label_data, mapping = aes(UMAP_1, UMAP_2, label = read_type), box.padding = 0.5,  force = 2)
    p2 <- p2 + facet_grid(~combined_cna_group)
    p2 <- p2 + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position = "none")
    file2write <- paste(dir_out_tmp2, snRNA_aliquot_id_tmp,".Malignant_Reclustered.CNA_", cna_gene1, "_", cna_gene2, ".Mutation_", mut_gene, ".", run_id, ".png", sep="")
    png(file2write, width = 3500, height = 800, res = 150)
    print(p2)
    dev.off()
  }
  
}

