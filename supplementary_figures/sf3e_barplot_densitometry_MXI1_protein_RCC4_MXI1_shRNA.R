# Yige Wu @WashU May 2022
## reference: https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/


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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
densitometry_df <- fread(data.table = F, input = "./Resources/Validation/Western_Blot/wb_densitometry_alltargets_normalized.10112022.v2.csv")

# set plot parameters -----------------------------------------------------
genes_plot <- c("MXI1")
lines_plot <- c("RCC4_scrambled", "RCC4_MXI1_C1")
linetexts_plot <- c("sh-NC", "sh-MXI1")
test_plot <- "t.test"
# test_plot <- "wilcox.test"

# preprocess --------------------------------------------------------------------
plotdata_df <- densitometry_df %>%
  mutate(Date = as.character(Date)) %>%
  mutate(Line = paste0(Parent_Line, "_", Construct_transfected)) %>%
  filter(Gene_symbol %in% genes_plot) %>%
  filter(Line %in% lines_plot) %>%
  filter(!is.na(Value_bytub_byscrambled))
plotdata_df <- data.frame(plotdata_df)
countsamples_bydate <- table(plotdata_df[, c("Date", "Line")]) %>% rowSums()
dates_filtered <- names(countsamples_bydate)[countsamples_bydate == length(lines_plot)]
plotdata_df <- plotdata_df %>%
  filter(Date %in% dates_filtered) %>%
  select(Line, Value_bytub_byscrambled)
# t.test(x = plotdata_df$Value_bytub_byscrambled[plotdata_df$Line == "RCC4_KLF9_C2"], y = rep(1, 4))
## p-value = 0.1051
## make colors
colors_byline <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[c(1:4, 6)]
names(colors_byline) <- c("RCC4_scrambled", "RCC4_KLF9_C2", "RCC4_MXI1_C1", "RCC4_MXI1_C2", "RCC4_KLF9_C3")
labels_x <- linetexts_plot; names(labels_x) <- lines_plot
## calculate limit for y-axis
plotdata_sum_df <- plotdata_df %>%
  group_by(Line) %>%
  summarise(y_plot = mean(Value_bytub_byscrambled),
            sd_plot = sd(Value_bytub_byscrambled))
ymax <- max(plotdata_sum_df$y_plot+plotdata_sum_df$sd_plot)
# test and plot -----------------------------------------------------------
stat.test <- compare_means(
  Value_bytub_byscrambled ~ Line, data = plotdata_df,
  method = test_plot, ref.group = "RCC4_scrambled"
)

p <- ggbarplot(data = plotdata_df, x = "Line", y = "Value_bytub_byscrambled", 
               add = c("mean_sd", "dotplot"), fill = "Line", color = "black", 
               error.plot = "upper_linerange")
p <- p + stat_pvalue_manual(stat.test, 
                            y.position = seq(ymax*1.05, ymax*(1+0.1*(length(lines_plot)-1)), length.out = (length(lines_plot)-1)), 
                            # label = "p = {p.format}",
                            label = "p = {signif(p, digits = 2)}",
                            label.size = 5)
p <- p + scale_fill_manual(values = colors_byline)
# p <- p + scale_x_discrete(labels = c("RCC4_scrambled" = "sh-NC", "RCC4_KLF9_C2" = "sh-KLF9"))
p <- p + scale_x_discrete(labels = labels_x)
# p <- p + ggtitle(label = paste0(genes_plot, " densitometry"))
p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, ymax*(1+0.1*(length(lines_plot))))) 
p <- p + theme_classic()
p <- p + ylab(label = paste0("Relative ", genes_plot, " level"))
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = "black", size = 12), 
               axis.text.y = element_text(color = "black", size = 12))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(color = "black", size = 12))
p
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".", paste0(lines_plot, collapse = "_"), ".", test_plot, ".densitometry_bytub_byscrambled.", "pdf")
# pdf(file2write, width = 1.4, height = 3.5, useDingbats = F) ## KLF9 - 2 lines
pdf(file2write, width = 1.75, height = 2.5, useDingbats = F) ## MXI1 - 3 lines
print(p)
dev.off()

## write source data
write.table(x = plotdata_df, file = paste0("~/Desktop/SF3e.SourceData.tsv"), quote = F, sep = "\t", row.names = F)
