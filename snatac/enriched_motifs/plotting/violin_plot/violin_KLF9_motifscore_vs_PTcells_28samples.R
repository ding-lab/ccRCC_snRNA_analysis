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
plotdata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/plotting/violin_plot/violin_selected_motifs_tumorcells_vs_PTcells_28samples/20220602.v1/MA1107.2.motif_score_by_cell.tsv", fill=TRUE)


# make colors -------------------------------------------------------------
color_tumorcell <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[1]

# plot --------------------------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(cell_group_text = ifelse(cell_group_plot %in% c("Tumor", "EMT tumor cells"), "Tumor cells", "PT cells"))
p <- ggviolin(data = plotdata_df, x = "cell_group_text", y = "motif_score", fill = "cell_group_text", color = NA,
               add = "boxplot", add.params = list(fill = "white", width = 0.15, color = "black"))
p <- p + scale_fill_manual(values = c("Tumor cells" = color_tumorcell, "PT cells" = color_pt))
p <- p + stat_compare_means(method = "t.test", label = "p.format", label.y = 6, label.x = 1.25)
p <- p + ylab(label = paste0("KLF9 ", " motif enrichment"))
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), 
               axis.title.y = element_text(size = 12, color = "black"), axis.text = element_text(color = "black", size = 12))
p
## write output
file2write <- paste0(dir_out, "KLF9", ".pdf")
pdf(file2write, width = 2.5, height = 2.5, useDingbats = F)
print(p)
dev.off()



