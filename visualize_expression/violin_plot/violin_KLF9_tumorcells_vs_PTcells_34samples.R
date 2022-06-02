# Yige Wu @WashU Jun 2022
## source activate signac

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "readxl",
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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data ------------------------------------------------------
plotdata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/visualize_expression/violin_plot/violin_selected_genes_tumorcells_vs_PTcells_34samples_katmai/20220602.v1/KLF9.expression_by_cell.tsv", fill=TRUE)

# plot --------------------------------------------------------------------
p <- ggviolin(data = plotdata_df, x = "cell_group_text", y = "exp_value", fill = "cell_group_plot", 
               add = "boxplot", add.params = list(fill = "white"))
p + stat_compare_means(method = "t.test", aes(label = paste0("p =", ..p.format..)))
## write output
file2write <- paste0(dir_out, gene_plot, ".pdf")
pdf(file2write, width = 2.5, height = 2.5, useDingbats = F)
print(p)
dev.off()



