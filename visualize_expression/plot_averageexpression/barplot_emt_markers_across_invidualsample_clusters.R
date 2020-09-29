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

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_byindividualcluster_bycellgroup7_byaliquot_on_katmai/20200917.v1/avgexp.SCT.bycellgroup.byaliquot.bycluster.31_aliquot_integration.20200917.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv", data.table = F)
## barcode 2 individual cluster id
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")


# specify the EMT samples ---------------------------------------------------
ids_aliquot_emt

# specify the gene to plot ------------------------------------------------
gene_plot <- "VIM"

# count cell number and filter clusters -----------------------------------
barcode2celltype_df <- merge(barcode2celltype_df, barcode2cluster_df, by.x = c("orig.ident", "individual_barcode"), by.y = c("aliquot", "individual_barcode"), all.x = T)
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_bycluster_bycellgroup_byaliquot = paste0(orig.ident, "_", seurat_cluster_id, "_",Cell_group7))
cellcount_bycluster_df <- barcode2celltype_df %>%
  select(id_bycluster_bycellgroup_byaliquot) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_bycluster_bycellgroup_byaliquot_original = ".") %>%
  mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = id_bycluster_bycellgroup_byaliquot_original, pattern = " |\\-", replacement = "."))

for (gene_plot in c("VIM", "LRP2")) {
  # format expression data --------------------------------------------------
  avgexp_df[1:5, 1:5]
  plot_data_df <- avgexp_df %>%
    filter(V1 == gene_plot) %>%
    melt() %>%
    mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
    mutate(aliquot = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,1]) %>%
    mutate(id_cluster = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,2]) %>%
    mutate(cellgroup = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,3])
  plot_data_df$Cell_count <- mapvalues(x = plot_data_df$id_bycluster_bycellgroup_byaliquot, from = cellcount_bycluster_df$id_bycluster_bycellgroup_byaliquot, to = as.vector(cellcount_bycluster_df$Freq))
  plot_data_df$Cell_count <- as.numeric(as.vector(plot_data_df$Cell_count))
  plot_data_df <- plot_data_df %>%
    dplyr::filter(Cell_count >= 30) %>%
    dplyr::filter(cellgroup %in% c("Tumor.cells", "Transitional.cells", "Tumor.like.cells")) %>%
    arrange(desc(value))
  plot_data_df$id_aliquot_wu <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
  ## make data to order the samples
  sample_ordered_df <- plot_data_df %>%
    group_by(id_aliquot_wu) %>%
    top_n(wt = value, n = 1)
  plot_data_df$id_aliquot_wu <- factor(x = plot_data_df$id_aliquot_wu, levels = sample_ordered_df$id_aliquot_wu)
  # write output ------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = id_aliquot_wu, y = value))
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p
  
  # write output ------------------------------------------------------------
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1200, height = 800, res = 150)
  print(p)
  dev.off()
}


