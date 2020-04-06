# Yige Wu @WashU March 2020
## make scatterplot to compare the xCell abundance values with snRNA cell group fractions

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
library(ggpubr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged xCell and snRNA cell group fractions
cellgroup_stat_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/deconvolution/merge_xCell_and_snRNA_cellgroup_fractions/20200320.v1/merged_xCell_value_and_snRNA_cellgroup_fractions.20200320.v1.tsv", data.table = F)

# loop for each cell group ------------------------------------------------
for (cellgroup.xcell_tmp in unique(cellgroup_stat_df$Cellgroup.xCell)) {
  ## make plot data
  plot_data_df <-  cellgroup_stat_df %>%
    filter(Cellgroup.xCell == cellgroup.xcell_tmp) %>%
    filter(!is.na(Value.xCell) & !is.na(CellgroupFrac.snRNA)) %>%
    mutate(x = CellgroupFrac.snRNA) %>%
    mutate(y = Value.xCell)
  
  for (method_tmp in c("pearson", "spearman")) {
    dir_out_tmp <- paste0(dir_out, method_tmp, "/")
    dir.create(dir_out_tmp)
    
    ## plot scatter plot
    p <- ggplot()
    p <- ggscatter(plot_data_df, x = "x", y = "y",
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
    )
    # Add correlation coefficient
    p <- p + stat_cor(method = method_tmp,  
                      label.x = quantile(x = plot_data_df$x, 0.5), 
                      label.y = quantile(x = plot_data_df$y, 0.95))
    p <- p + xlab("Cell Group Fraction Estimated from snRNA Data")
    p <- p + ylab("Cell Group Abundance Estimated from Bulk RNA Data (xCell)")
    p <- p + ggtitle(label = paste0(cellgroup.xcell_tmp, " Fraction Comparison"), subtitle = "snRNA-based Cell Type Fraction vs. xCell")
    p
    
    ## save scatterplot
    plot_path <- paste0(dir_out_tmp, "scatterplot_", cellgroup.xcell_tmp, "_from_bulk_and_snRNA.", run_id, ".png")
    png(filename = plot_path, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
  }
}