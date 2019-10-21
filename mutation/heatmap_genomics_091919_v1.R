# Yige Wu @WashU Sep 2019
## for plotting the genomics landscape of integrated single cell data

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
sample_ids <- c("CPT0086820004", "CPT0075130004")
run_id <- paste0(sample_ids, collapse = "_")

# input 10Xmapping results ------------------------------------------------
sup_tab_10Xmapping <- NULL
for (sample_id_tmp in sample_ids) {
  dir_10Xmapping_tmp <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/", 
                               sample_id_tmp, "/", sample_id_tmp, "/10Xmapping/")
  file_names_tmp <- list.files(dir_10Xmapping_tmp)
  file_names_tmp <- file_names_tmp[grepl(pattern = "mapping_heatmap", x = file_names_tmp)]
  file_names_tmp
  for (file_name_tmp in file_names_tmp) {
    path_file_tmp <- paste0(dir_10Xmapping_tmp, file_name_tmp)
    file_tmp <- fread(input = path_file_tmp, data.table = F)
    file_tmp <- file_tmp %>%
      mutate(Mutation = Mutatation) %>%
      mutate(gene_symbol = str_split_fixed(string = Mutation, pattern = "-", n = 3)[,1]) %>%
      filter(gene_symbol %in% SMGs[["CCRCC"]])
    file_tmp$Mutatation <- NULL
    file_tmp.m <- melt(data = file_tmp, id.vars = c("Mutation", "gene_symbol"))
    file_tmp.m <- file_tmp.m %>%
      mutate(cluster_name = str_split_fixed(string = file_name_tmp, pattern = "mapping_heatmap_|.txt", n = 3)[,2]) %>%
      mutate(allele_type = str_split_fixed(string = Mutation, pattern = "\\-", n = 3)[,3]) %>%
      mutate(aliquot = sample_id_tmp) %>%
      filter(!is.na(value))
    file_tmp.m
    sup_tab_10Xmapping <- rbind(file_tmp.m, sup_tab_10Xmapping)
  }
}

# plot heatmap ------------------------------------------------------------
tab2p <- sup_tab_10Xmapping %>%
  filter(allele_type == "Var") %>%
  dplyr::group_by(aliquot, Mutation) %>%
  summarize(num_reads = sum(value))

mat2p <- dcast(data = tab2p, formula = Mutation ~ aliquot)
rownames(mat2p) <- mat2p$Mutation
mat2p$Mutation <- NULL
mat2p <- as.matrix(x = mat2p)
mat2p

## get heatmap colors
col_fun = colorRamp2(c(0, 8), c("white", "red"))
col_fun(seq(0, 8))

p <- Heatmap(mat2p, 
             cluster_rows = F,
             cluster_columns = F,
             col = col_fun,
             name = "Num.Var.Reads",
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (!is.na(mat2p[i, j]))
               grid.text(sprintf("%d", mat2p[i, j]), x, y, gp = gpar(fontsize = 10))
             },
             na_col = "white")
version_tmp <- 1
pdf(file = paste0(makeOutDir(), "10Xmapping_only.", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".pdf"), width = 4, height = 5)
print(p)
dev.off()