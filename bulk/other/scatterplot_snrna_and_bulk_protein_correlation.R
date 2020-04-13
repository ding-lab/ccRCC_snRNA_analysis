# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggpubr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input denpendencies -----------------------------------------------------
spearman_result_df <- fread(input = "./Resources/Analysis_Results/bulk/other/correlate_snrna_and_bulk_protein/20200411.v1/spearman_result.20200411.v1.tsv", data.table = F)
## input differential gene result
deg_roc_df <- fread(input = "./Resources/Analysis_Results/findmarkers/filter_celltypemarkers_for_bulk_correlation/20200411.v1/celltypemarkers_filtered_for_bulk_correlation.20200411.v1.tsv", data.table = F)
## input snRNA averaged by cell type
avgexp_bycelltype_df <- fread(input = "./Resources/Analysis_Results/average_expression/adjust_averageexpression_for_bulk_correlation/20200411.v1/avgexp_bycelltype.long_data_frame.20200411.v1.tsv", data.table = F)
## load bulk protein data
protein_df <- fread("./Resources/Bulk_Processed_Data/Protein/CCRCC_PRO_tumor_PGDAC_MD_MAD_partID.txt", data.table = F)
## load id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input barcode 2 cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv", data.table = F)


# count number of cells per celltype per aliquot --------------------------
count_bycelltype_byaliquot_df <- data.frame(table(barcode2celltype_df[, c("Cell_type.shorter", "orig.ident")]))
unique(count_bycelltype_byaliquot_df$Cell_type.shorter)
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$Cell_type.shorter, pattern = '[/ +]', replacement = ".")
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$celltype, pattern = '\\-', replacement = ".")

# filter and map case id --------------------------------------------------
## only keep the original piece
avgexp_bycelltype_df <- avgexp_bycelltype_df %>%
  filter(aliquot %in% idmetadata_df$Aliquot.snRNA[idmetadata_df$Is_discovery_set == T & idmetadata_df$Sample_Type == "Tumor"])
## map case id
avgexp_bycelltype_df$id_case <- mapvalues(x = avgexp_bycelltype_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

# filter to key genes -----------------------------------------------------
genes_celltype <- spearman_result_df$gene_symbol[spearman_result_df$pvalue < 0.05 & spearman_result_df$rho > 0 & spearman_result_df$celltype != "All_Cells"]
genes_allcells<- spearman_result_df$gene_symbol[spearman_result_df$pvalue < 0.05 & spearman_result_df$rho > 0 & spearman_result_df$celltype == "All_Cells"]
genes_plot <- intersect(genes_celltype, genes_allcells)
genes_plot

# plot --------------------------------------------------------------------
for (gene_tmp in genes_plot) {
  celltypes_plot <- unique(spearman_result_df$celltype[spearman_result_df$gene_symbol == gene_tmp & spearman_result_df$pvalue < 0.05])
  exp_df <- avgexp_bycelltype_df %>%
    filter(gene_symbol == gene_tmp) %>%
    rename(snrna = avg_exp)
  ## use the case id to extract the protein
  idcase2test <- exp_df$id_case
  protein_tmp <- protein_df[protein_df$Gene == gene_tmp, idcase2test]
  protein_tmp <- as.numeric(unlist(protein_tmp))
  ## add protein
  exp_df$protein <- protein_tmp
  ## filter out exp values calculated from too few cells
  exp_df <- merge(exp_df, count_bycelltype_byaliquot_df, by.x = c("aliquot", "celltype"), by.y = c("orig.ident", "celltype"), all.x = T)
  exp_df <- exp_df %>%
    filter((celltype == "All_Cells") | (Freq >= 10))
  
  for (celltype_tmp in celltypes_plot) {
    ## make plot for showing cell type expression compared to other cell types
    plot_data_df <- exp_df %>%
      filter(!(celltype %in% c("All_Cells", "Unknown"))) %>%
      mutate(x = protein) %>%
      mutate(y = snrna) %>%
      mutate(color = ifelse(celltype %in% celltype_tmp, celltype_tmp, "Other_Cell_Types"))
    colors_bycelltype <- c("red", "grey")
    names(colors_bycelltype) <- c(celltype_tmp, "Other_Cell_Types")
    p <- ggplot()
    p <- p + geom_point(data = plot_data_df, mapping = aes(x = x, y = y, color = color), alpha = 0.6)
    p <- p + scale_color_manual(values = colors_bycelltype)
    p <- p + xlab(label = "Bulk Protein Level (Relative to Reference)")
    p <- p + ylab(label = "snRNA Expression Level (Averaged by Cell Type)")
    p <- p + theme(legend.position = "top")
    p
    file2write <- paste0(dir_out, gene_tmp, ".snRNA_protein.", celltype_tmp, "_highlighted.", run_id, ".png")
    png(filename = file2write, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
    ## make scatterplot for showing correlation
    plot_data_df <- exp_df %>%
      filter(celltype == celltype_tmp) %>%
      mutate(x = protein) %>%
      mutate(y = snrna)
    p <- ggplot()
    p <- ggscatter(data = plot_data_df, x = "x", y = "y",
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
    )
    p <- p + stat_cor(method = "spearman", label.x = 0, label.y = quantile(x = plot_data_df$y, probs = 0.99))
    p <- p + xlab(label = "Bulk Protein Level (Relative to Reference)")
    p <- p + ylab(label = "snRNA Expression Level (Averaged by Cell Type)")
    p <- p + ggtitle(label = paste0(gene_tmp, " snRNA Expression in ", celltype_tmp, "~ Bulk Protein Level"))
    p
    file2write <- paste0(dir_out, gene_tmp, ".snRNA_protein.", celltype_tmp, "_scatterplot.", run_id, ".png")
    png(filename = file2write, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
  }
}

