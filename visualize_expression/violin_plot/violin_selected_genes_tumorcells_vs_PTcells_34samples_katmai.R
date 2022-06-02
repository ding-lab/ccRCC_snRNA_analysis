# Yige Wu @WashU Jun 2022
## source activate signac

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
  "ggplot2",
  "ggpubr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input data ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20211005.v1//ccRCC.34samples.Merged.20211005.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")

# set parameters ----------------------------------------------------------
gene_plot <- "KLF9"

# get plot data -----------------------------------------------------------
exp_mat <- srat[["SCT"]]@data
exp_df <- data.frame(barcode_integrated = colnames(exp_mat), exp_value = as.vector(exp_mat[gene_plot,]))
exp_df$orig.ident <- mapvalues(x = exp_df$barcode_integrated, from = rownames(srat@meta.data), to = as.vector(srat@meta.data$orig.ident))
rm(exp_mat)
rm(srat)
exp_df <- exp_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode_integrated, pattern = "_", n = , 2)[,1]) %>%
  mutate(barcode_uniq = paste0(orig.ident, "_", barcode_individual))
## add cell type annotation
exp_df <- merge(x = exp_df, 
                y = barcode2celltype_df %>%
                  mutate(barcode_uniq = paste0(orig.ident, "_", individual_barcode)) %>%
                  mutate(cell_group_plot = Cell_group_w_epithelialcelltypes) %>%
                  select(barcode_uniq, cell_group_plot), by = c("barcode_uniq"), all.x = T)
## filter
plotdata_df <- exp_df %>%
  filter(cell_group_plot %in% c("Tumor cells", "Proximal tubule"))
plotdata_df$cell_group_text <- mapvalues(x = plotdata_df$cell_group_plot, from = c("Tumor cells", "Proximal tubule"), to = c("Tumor\ncells", "PT\ncells"))
# plotdata_df$cell_group_plot <- factor(x = plotdata_df$cell_group_plot, levels = c("Tumor cells", "Proximal tubule"))
plotdata_df$cell_group_text <- factor(x = plotdata_df$cell_group_text, levels = c("Tumor\ncells", "PT\ncells"))

## write output
file2write <- paste0(dir_out, gene_plot, ".expression_by_cell.tsv")
write.table(x = plotdata_df, file = file2write, quote = F, sep = "\t", row.names = F)

# plot --------------------------------------------------------------------
p <- ggviolin(data = plotdata_df, x = "cell_group_text", y = "exp_value", fill = "cell_group_plot", 
               add = "boxplot", add.params = list(fill = "white"))
p + stat_compare_means(method = "t.test", label = aes(label = paste0("p =", ..p.format..)))
## write output
file2write <- paste0(dir_out, gene_plot, ".pdf")
pdf(file2write, width = 2.5, height = 2.5, useDingbats = F)
print(p)
dev.off()



