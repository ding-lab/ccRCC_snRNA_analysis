# Yige Wu @WashU Feb 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution is too confusing

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
### set output subdirectories
dir_out1 <- paste0(dir_out, "All_CNVs/")
dir.create(dir_out1)
### set output subdirectories
dir_out2 <- paste0(dir_out, "Expected_CNVs/")
dir.create(dir_out2)
### set output subdirectories
dir_out3 <- paste0(dir_out, "Expected_CNVs_Higher_Than_NAT/")
dir.create(dir_out3)

# input dependencies ------------------------------------------------------
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200505.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200505.v1.tsv", data.table = F)
table(cnv_3state_count_aliquots$tumor_subcluster)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# process known cnv genes --------------------------------------------------
knowncnvgenes_df <- knowncnvgenes_df %>%
  mutate(Arm = ifelse(grepl(pattern = "p", x = Cytoband), "p", "q")) %>%
  mutate(Chr = str_split_fixed(string = Cytoband, pattern = Arm, n = 2)[,1])

# make matrix for heatmap body --------------------------------------------
heatmapbody_long_df <- cnv_3state_count_aliquots
heatmapbody_long_df$gene_expected_state <- mapvalues(x = heatmapbody_long_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
## add aliquot.wu
heatmapbody_long_df$aliquot.wu <- mapvalues(x = heatmapbody_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## filter out copy neutral ones
heatmapbody_long_df <- heatmapbody_long_df %>%
  filter(cna_3state != "Neutral") %>%
  filter(gene_expected_state == cna_3state) %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "C", tumor_subcluster))
## make matrix
heatmapbody_wide_df <- dcast(data = heatmapbody_long_df, formula = gene_symbol ~ id_aliquot_cluster, value.var = "Fraction")
heatmapbody_mat <- as.matrix(heatmapbody_wide_df %>%
                               select(-gene_symbol))
rownames_vec <- heatmapbody_wide_df$gene_symbol
rownames(heatmapbody_mat) <- rownames_vec
heatmapbody_mat[1:5, 1:5]

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = rownames_vec, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
row_split_vec
## order cytoband
knowncnvgenes_df <- knowncnvgenes_df %>%
  mutate(Chr = factor(x = Chr, levels = 1:22)) %>%
  arrange(Chr, Arm, Cytoband)
row_split_factor <- factor(x = row_split_vec, levels = unique(knowncnvgenes_df$Cytoband))

# make heatmap ------------------------------------------------------------
Heatmap(heatmapbody_mat, 
        name = "Fraction", 
        # col = col_fun, 
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width*heatmapbody_mat[i, j], height = height*heatmapbody_mat[i, j], 
                    gp = gpar(col = "grey", fill = NA))
          # grid.circle(x = x, y = y, r = abs(heatmapbody_mat[i, j])/2 * min(unit.c(width, height)), 
          #             gp = gpar(fill = "red", col = NA))
        }, 
        row_split = row_split_vec,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_names = T, 
        show_column_names = T)
