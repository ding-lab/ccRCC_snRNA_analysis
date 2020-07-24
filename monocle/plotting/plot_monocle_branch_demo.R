library(monocle)
library(tidyverse)
library(pheatmap)
my_branchplot <- function(cds_subset, 
                          branch_point = 1, 
                          branch_states = NULL, 
                          branch_labels = c("Cell fate 1", "Cell fate 2"), 
                          cluster_rows = TRUE, 
                          hclust_method = "ward.D2", 
                          num_clusters = 6, 
                          hmcols = NULL, 
                          branch_colors = c("#979797", "#F05662", "#7990C8"), 
                          add_annotation_row = NULL, 
                          add_annotation_col = NULL, 
                          show_rownames = FALSE, 
                          use_gene_short_name = TRUE, 
                          scale_max = 3, 
                          scale_min = -3, 
                          norm_method = c("log", "vstExprs"), 
                          trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", 
                          return_heatmap = FALSE, 
                          cores = 1, ...) 
{
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
  
  
  
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs") {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
  }
  else if (norm_method == "log") {
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
  df <- df %>% 
    select(Time,State,Branch,Pseudotime)
  
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
  pheatmap(heatmap_matrix, 
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
           filename = NA, 
           #breaks = bks, 
           color = hmcols)
}
### test it if it worked
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))

lung <- load_lung()
my_branchplot(cds_subset = lung,
              branch_point = 1,
              hmcols = mypalette,
              show_rownames = T)
