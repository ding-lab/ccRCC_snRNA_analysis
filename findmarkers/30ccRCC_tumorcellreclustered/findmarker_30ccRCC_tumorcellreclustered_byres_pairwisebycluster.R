# Yige Wu @WashU Mar 2022
## https://satijalab.org/seurat/archive/v3.0/integration.html
## also used references

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "future",
  "future.apply",
  "ggplot2"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out_parent <- makeOutDir_katmai(path_this_script)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
logfc.threshold.run <- 0.25
min.pct.run <- 0.1
min.diff.pct.run <- 0
## input the barcode-to-cluster results
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")

# process -----------------------------------------------------------------
resolutions_process <- colnames(barcode2cluster_df)[grepl(pattern = "integrated_snn_res", x = colnames(barcode2cluster_df))]
resolutions_process <- gsub(x = resolutions_process, pattern = "integrated_snn_res\\.", replacement = "")
results_sup_df <- NULL
for (resolution_tmp in c("1", "2")) {
# for (resolution_tmp in resolutions_process) {
  path_markers <- paste0(dir_out_parent, "res.", resolution_tmp, "tumorcellsreclustered.pairwisebycluster.markers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
  if (file.exists(path_markers)) {
    results_df <- fread(data.table = F, input = path_markers)
    print("Markers for resolution ", resolution_tmp, "exists, reading!\n")
  } else {
    print("Markers for resolution ", resolution_tmp, "doesn't exist, running FindMarkers!\n")
    
    srat@meta.data$cluster_test <- mapvalues(x = rownames(srat@meta.data), from = barcode2cluster_df$barcode, to = as.vector(barcode2cluster_df[, paste0("integrated_snn_res.", resolution_tmp)]))
    Idents(srat) <- "cluster_test"
    clusters <- unique(Idents(srat))
    pairwise <- combn(clusters, 2)
    
    results_df <- NULL
    for (i in 1:ncol(pairwise)) {
      markers <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = pairwise[1, i], ident.2 = pairwise[2, i], only.pos = T,
                             min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
      markers$gene_symbol <- rownames(markers)
      markers$ident.1 <- pairwise[1, i]
      markers$ident.2 <- pairwise[2, i]
      results_df <- rbind(results_df, markers)
    }
    write.table(x = results_df, file = path_markers, quote = F, sep = "\t", row.names = F)
  }
  results_sup_df <- rbind(results_sup_df, results_df)
  
  # plot --------------------------------------------------------------------
  plotdata_df <- results_df %>%
    filter(p_val_adj < 0.05) %>%
    mutate(diff_pct = (pct.1 - pct.2)) %>%
    # filter(diff_pct >= 0.1) %>%
    filter(diff_pct >= 0) %>%
    mutate(x_plot = ifelse(ident.1 > ident.2, ident.1, ident.2)) %>%
    mutate(y_plot = ifelse(ident.1 > ident.2, ident.2, ident.1)) %>%
    group_by(x_plot, y_plot) %>%
    summarise(number_degs = n())
  plotdata_df$x_plot <- factor(plotdata_df$x_plot)
  plotdata_df$y_plot <- factor(plotdata_df$y_plot)
  
  p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
  p <- p + geom_tile(mapping = aes(fill = number_degs_brks))
  p <- p + geom_text(mapping = aes(label = number_degs))
  p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                 midpoint = 20, space = "Lab",
                                 name="Number of\nDEGs")
  p <- p + xlab("Cluster #1") + ylab("Cluster #2")
  p <- p + theme_minimal()
  
  file2write <- paste0(dir_out, "Res", resolution_tmp, ".Number_of_DEGs.png")
  png(file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}


# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumorcellsreclustered.pairwisebycluster.markers.logfcthreshold.", 
                     logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run,
                     ".byresolution", ".tsv")
write.table(x = results_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)


