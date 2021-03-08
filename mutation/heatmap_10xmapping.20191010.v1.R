# Yige Wu @WashU Oct 2019
## for plotting the mutation landscape of integrated single cell data

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# specify aliquot ids to process --------------------------------------------
aliquot_ids <- c("CPT0075130004_notFACS", "CPT0086820004_notFACS", "CPT0075140002", "CPT0001260013", "CPT0086350004")

# input 10Xmapping results ------------------------------------------------
sup_tab_10Xmapping <- NULL
for (aliquot_id_tmp in aliquot_ids) {
  dir_10Xmapping_tmp <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/10Xmapping/outputs/", 
                               aliquot_id_tmp, "/")
  file_names_tmp <- list.files(dir_10Xmapping_tmp)
  file_names_tmp <- file_names_tmp[grepl(pattern = "mapping_heatmap", x = file_names_tmp)]
  file_names_tmp
  for (file_name_tmp in file_names_tmp) {
    path_file_tmp <- paste0(dir_10Xmapping_tmp, file_name_tmp)
    file_tmp <- fread(input = path_file_tmp, data.table = F)
    file_tmp <- file_tmp %>%
      mutate(Mutation = Mutatation) %>%
      mutate(gene_symbol = str_split_fixed(string = Mutation, pattern = "-", n = 3)[,1]) %>%
      dplyr::filter(gene_symbol %in% SMGs[["CCRCC"]])
    file_tmp$Mutatation <- NULL
    file_tmp.m <- melt(data = file_tmp, id.vars = c("Mutation", "gene_symbol"))
    file_tmp.m <- file_tmp.m %>%
      mutate(cluster_name = str_split_fixed(string = file_name_tmp, pattern = "mapping_heatmap_|.txt", n = 3)[,2]) %>%
      mutate(allele_type = str_split_fixed(string = Mutation, pattern = "\\-", n = 3)[,3]) %>%
      mutate(aliquot = aliquot_id_tmp) %>%
      dplyr::filter(!is.na(value))
    file_tmp.m
    sup_tab_10Xmapping <- rbind(file_tmp.m, sup_tab_10Xmapping)
  }
}

# plot heatmap ------------------------------------------------------------
tab2p <- sup_tab_10Xmapping %>%
  mutate(mutation_text = str_split_fixed(string = Mutation, pattern = "-Var", n = 2)[,1]) %>%
  dplyr::filter(allele_type == "Var") %>%
  dplyr::group_by(aliquot, mutation_text) %>%
  dplyr::summarize(num_reads = sum(value))

mat2p <- dcast(data = tab2p, formula = aliquot ~ mutation_text, value.var = "num_reads")
rownames(mat2p) <- mat2p$aliquot
mat2p$aliquot <- NULL
mat2p <- as.matrix(x = mat2p)
mat2p

ref_count_mat <- dcast(data = sup_tab_10Xmapping %>%
                         mutate(mutation_text = str_split_fixed(string = Mutation, pattern = "-Ref", n = 2)[,1]) %>%
                         dplyr::filter(allele_type == "Ref") %>%
                         dplyr::group_by(aliquot, mutation_text) %>%
                         dplyr::summarize(num_reads = sum(value)), 
                       formula = aliquot ~ mutation_text, value.var = "num_reads")
ref_count_mat
rownames(ref_count_mat) <- ref_count_mat$aliquot
ref_count_mat$aliquot <- NULL
ref_count_mat <- as.matrix(x = ref_count_mat)
ref_count_mat

## get heatmap colors
col_fun = colorRamp2(c(0, 8), c("white", "red"))
col_fun(seq(0, 8))

## get colume text

p <- Heatmap(mat2p, 
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F, 
             top_annotation = HeatmapAnnotation(foo = anno_text(colnames(mat2p), 
                                                                rot = 90,
                                                                # location = unit(1, 'npc'),
                                                                gp = gpar(fontsize = 12))),
             col = col_fun,
             name = "No.Var/Ref",
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (!is.na(mat2p[i, j]))
               grid.text(sprintf("%s", paste0(mat2p[i, j], "/", ref_count_mat[i, j])), x, y, gp = gpar(fontsize = 10))
             },
             na_col = "white")

png(file = paste0(dir_out, "10Xmapping", ".", run_id, ".png"), width = 1200, height = 500, res = 150)
print(p)
dev.off()
