# Yige Wu @WashU May 2021

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
## input motif correlation
sharedgene_frac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/calculate_motif_shared_degs/20210513.v1/SharedGeneFraction.tsv")


# make plot data ----------------------------------------------------------
plotdata_mat <- as.matrix(sharedgene_frac_df[,-1])
rownames(plotdata_mat) <- sharedgene_frac_df$TF_name


# make row annotation -----------------------------------------------------
highlight_motif <- "RBPJ"
column_order <- sharedgene_frac_df[,highlight_motif]
names(column_order) <- sharedgene_frac_df$TF_name
top_values <- head(column_order[order(column_order, decreasing = T)], n = 10)
ha = rowAnnotation(foo = anno_mark(at = which(column_order %in% top_values), labels = names(column_order)[column_order %in% top_values], labels_gp = gpar(fontsize = 40)))

# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plotdata_mat,
             # col = col_fun, 
             ## row
             # row_split = factor_case_ids, cluster_row_slices = F,
             # row_title = NULL,
             # # row_title_rot = 0, row_title_gp = gpar(fontsize = 80),
             show_row_dend = F,
             # row_gap = unit(0, "mm"),
             right_annotation = ha,
             show_row_names = F,
             # 
             ## column
             # column_split = factor_case_ids, cluster_column_slices = F, 
             # column_title = NULL,
             # # column_title_side = "bottom", column_title_rot = 90, column_title_gp = gpar(fontsize = 80),
             show_column_dend = F,
             # column_gap = unit(0, "mm"),
             # # top_annotation = col_anno_obj2,
             # top_annotation= col_anno_obj1,
             # show_column_names = F,
             show_heatmap_legend = T)
file2write <- paste0(dir_out, "Pariwise_Correlation", ".pdf")
pdf(file2write,
    width = 60, height = 60, useDingbats = F)
draw(p)
dev.off()
