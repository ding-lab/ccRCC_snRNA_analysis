# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
library(pheatmap)
library(tidyverse)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input monocle object
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3N-01200/combined_subset_pseudotime_qval_1e-10.rds")
## input deg table
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/monocle/deg/unite_pseudotime_degs_across_3samples/20200728.v1/PT_and_TumorCells_ByCellType.Pseudotim_DE_genes.Combine3Samples.20200728.v1.tsv")

# prepare parameters --------------------------------------------------
deg_filtered <- deg_df %>%
  filter(use_for_ordering) %>%
  filter(qval<0.1)
deg_count_df <- deg_filtered %>%
  select(gene_short_name) %>%
  table() %>%
  data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(gene_symbol = ".")
## specify genes to plot
genes_plot <- as.vector(deg_count_df$gene_symbol[deg_count_df$Freq == 3])
# genes_plot <- genes_plot[1:20]
## make color palette
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))

# specify parameters ------------------------------------------------------
## subset the monocle object
cds_subset <- obj_monocle[genes_plot,]
branch_point = 1
branch_states <- NULL
cores = 1
trend_formula = "~sm.ns(Pseudotime, df=3) * Branch"
scale_max = 3
scale_min = -3
hmcols <- mypalette
norm_method <- "vstExprs"
hclust_method = "ward.D2"
show_rownames <- T
num_clusters = 4
# get data matrix for the heatmap ------------------------------------------------------
pData(cds_subset)$Pseudotime <- (pData(cds_subset)$Pseudotime*100)/max(pData(cds_subset)$Pseudotime) 
new_cds <- buildBranchCellDataSet(cds_subset, 
                                  branch_states = branch_states, 
                                  branch_point = branch_point, 
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo

if (is.null(branch_states)) {
  progenitor_state <- subset(pData(cds_subset), Pseudotime == 
                               0)[, "State"]
  branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
}
new_cds.pdata <- pData(new_cds)
new_cds.pdata$Cell_ID <- rownames(new_cds.pdata)
col_gap_ind <- table(new_cds.pdata$Branch)[[1]]

pseudotime.A <- new_cds.pdata %>%
  filter(Branch==unique(as.character(pData(new_cds)$Branch))[1]) %>%
  select(Cell_ID,Pseudotime)

pseudotime.A <- pseudotime.A[order(pseudotime.A$Pseudotime),]

pseudotime.B <- new_cds.pdata %>%
  filter(Branch==unique(as.character(pData(new_cds)$Branch))[2]) %>%
  select(Cell_ID,Pseudotime) 

pseudotime.B <- pseudotime.B[order(pseudotime.B$Pseudotime),]

newdataA <- data.frame(Pseudotime = pseudotime.A$Pseudotime, 
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))

newdataB <- data.frame(Pseudotime = pseudotime.B$Pseudotime, 
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))

BranchAB_exprs <- genSmoothCurves(new_cds[, ], 
                                  cores = cores, 
                                  trend_formula = trend_formula, 
                                  relative_expr = T, 
                                  new_data = rbind(newdataA,newdataB))

colnames(BranchAB_exprs) <- c(pseudotime.A$Cell_ID,pseudotime.B$Cell_ID)
BranchA_exprs <- BranchAB_exprs[, 1:col_gap_ind]
BranchB_exprs <- BranchAB_exprs[, (col_gap_ind+1):ncol(BranchAB_exprs)]

# norm_method <- match.arg(norm_method)
if (norm_method == "vstExprs") {
  BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
  BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
} else if (norm_method == "log") {
  BranchA_exprs <- log10(BranchA_exprs + 1)
  BranchB_exprs <- log10(BranchB_exprs + 1)
}
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], 
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, 
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), 
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == 
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
# exp_rng <- range(heatmap_matrix)
# bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
# if (is.null(hmcols)) {
#   hmcols <- blue2green2red(length(bks) - 1)
# }

df <- new_cds.pdata[-col_gap_ind,]
### for selecting which features go into the the column annotation
df <- df %>% 
  select(Aliquot.snRNA.WU, State, Branch, Cell_type.detailed, Pseudotime)

###only clustering not plot
ph <- pheatmap(heatmap_matrix, 
               useRaster = T, 
               cluster_cols = FALSE, 
               cluster_rows = TRUE, 
               show_rownames = F, 
               show_colnames = F, 
               annotation_col = df,
               gaps_col = col_gap_ind,
               clustering_distance_rows = row_dist, 
               clustering_method = hclust_method, 
               cutree_rows = num_clusters, 
               silent = F, ###not plot
               filename = NA, 
               #breaks = bks, 
               color = hmcols)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))

# plot --------------------------------------------------------------------
p <- pheatmap(heatmap_matrix, 
         useRaster = T, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         show_rownames = show_rownames, 
         show_colnames = F, 
         annotation_col = df,
         annotation_row = annotation_row,
         gaps_col = col_gap_ind,
         clustering_distance_rows = row_dist, 
         clustering_method = hclust_method, 
         cutree_rows = num_clusters, 
         silent = F, ###not plot
         #breaks = bks, 
         color = hmcols)


# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "Pseudotie_DEG.png")
save_pheatmap_png(x = p, filename = file2write, width = 2000, height = 2200, res = 150)
file2write <- paste0(dir_out, "Pseudotie_DEG.pdf")
save_pheatmap_pdf(x = p, filename = file2write, width = 10, height = 15)


