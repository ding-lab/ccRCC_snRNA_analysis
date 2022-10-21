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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
# densitometry_df <- readxl::read_xlsx(path = "./Validation/Western_Blot/wb densitometry_05162022_pl.xlsx", sheet = "mxi1_klf9_normalized_05162022")
# densitometry_df <- readxl::read_xlsx(path = "./Validation/Western_Blot/wb densitometry_05162022_pl.xlsx", sheet = "mxi1_klf9_normalized_v2")
densitometry_df <- readxl::read_xlsx(path = "./Validation/Western_Blot/wb densitometry_05172022_pl_v2.xlsx", sheet = "alltargets_normalized")
densitometry_df <- fread(data.table = F, input = "./Validation/Western_Blot/wb_densitometry_alltargets_normalized.csv")

# set plot parameters -----------------------------------------------------
genes_plot <- c("CP")
genes_plot <- c("KLF9")
# genes_plot <- c("KLF9")
# genes_plot <- c("HK2")
# genes_plot <- c("PFKP")
# genes_plot <- c("PKM2")
# lines_plot <- c("RCC4_scrambled", "RCC4_KLF9_C2", "RCC4_KLF9_C3")
lines_plot <- c("RCC4_scrambled", "RCC4_KLF9_C2")
# genes_plot <- c("MXI1")
# genes_plot <- c("CP")
# genes_plot <- c("HK2")
# lines_plot <- c("RCC4_scrambled", "RCC4_MXI1_C1", "RCC4_MXI1_C2")
# lines_plot <- c("RCC4_scrambled", "RCC4_MXI1_C2")
test_plot <- "t.test"
# test_plot <- "wilcox.test"

# preprocess --------------------------------------------------------------------
plotdata_df <- densitometry_df %>%
  mutate(Date = as.character(Date)) %>%
  mutate(Line = paste0(Parent_Line, "_", Construct_transfected)) %>%
  filter(Gene_symbol %in% genes_plot) %>%
  filter(Line %in% lines_plot) %>%
  filter(!is.na(Value_bytub_byscrambled))
countsamples_bydate <- table(plotdata_df %>%
        select(Date, Line)) %>% rowSums()
dates_filtered <- names(countsamples_bydate)[countsamples_bydate == length(lines_plot)]
plotdata_df <- plotdata_df %>%
  # filter(Date != "2022-04-23") %>%
  # filter(Date != "2022-05-12") %>%
  filter(Date %in% dates_filtered)
t.test(x = plotdata_df$Value_bytub_byscrambled[plotdata_df$Line == "RCC4_KLF9_C2"], y = rep(1, 4))
## p-value = 0.1051
## make colors
colors_byline <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[c(1:4, 6)]
names(colors_byline) <- c("RCC4_scrambled", "RCC4_KLF9_C2", "RCC4_MXI1_C1", "RCC4_MXI1_C2", "RCC4_KLF9_C3")
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
               add = "mean_sd", fill = "Line", color = "black", 
               error.plot = "upper_linerange")
p <- p + stat_pvalue_manual(stat.test, 
                            y.position = seq(ymax*1.05, ymax*(1+0.1*(length(lines_plot)-1)), length.out = (length(lines_plot)-1)), 
                            # label = "p = {p.format}",
                            label = "p = {signif(p, digits = 2)}",
                            label.size = 5)
p <- p + scale_fill_manual(values = colors_byline)
p <- p + scale_x_discrete(labels = c("RCC4_scrambled" = "sh-NC", "RCC4_KLF9_C2" = "sh-KLF9"))
# p <- p + ggtitle(label = paste0(genes_plot, " densitometry"))
# p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, 1.15)) ## KLF9
# p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, 1.3)) ## MXI1 expression for MXI1 lines - 3 lines
p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, ymax*(1+0.1*(length(lines_plot))))) ## CP expression for MXI1 lines - 3 lines
# p <- p  + scale_y_continuous(expand = c(0,0), limits = c(0, ymax*(1+0.2*(length(lines_plot))))) ## CP expression for MXI1 lines - 3 lines
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
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".", paste0(lines_plot, collapse = "_"), ".", test_plot, ".densitometry_bytub_byscrambled.", "png")
png(file2write, width = 200, height = 500, res = 150)
print(p)
dev.off()

# stat.test <- compare_means(
#   Value_bytub ~ Line, data = plotdata_df,
#   method = test_plot, ref.group = "RCC4_scrambled"
# )
# 
# p <- ggbarplot(data = plotdata_df, x = "Line", y = "Value_bytub", add = "mean_sd", fill = "Construct_transfected")
# p <- p + stat_pvalue_manual(stat.test, 
#                             y.position = seq(1.05, 1.5, 0.05)[1:(length(lines_plot)-1)],
#                             label = "p.signif")
# # p <- p + scale_fill_manual(values = c("scrambled" = "grey", ""))
# p <- p + ggtitle(label = paste0(genes_plot, " densitometry"))
# # p <- p + ylim(c(0, 1.25))
# p <- p + theme_classic()
# p <- p + ylab(label = "% to control")
# p <- p + theme(legend.position = "none")
# p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
# p
# file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".", paste0(lines_plot, collapse = "_"), ".", test_plot, ".densitometry_bytub.", "pdf")
# pdf(file2write, width = 3, height = 3.5, useDingbats = F)
# print(p)
# dev.off()
# file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".", paste0(lines_plot, collapse = "_"), ".", test_plot, ".densitometry_bytub.", "png")
# png(file2write, width = 400, height = 500, res = 150)
# print(p)
# dev.off()


# backup ------------------------------------------------------------------
# stat.test <- compare_means(
#   Value_bytub_byscrambled ~ Line, data = plotdata_df,
#   method = test_plot,
# )

# plotdata_sum_df <- plotdata_df %>%
#   group_by(Line) %>%
#   summarise(y_plot = mean(Value_bytub_byscrambled),
#             sd_plot = sd(Value_bytub_byscrambled))
# plotdata_sum_df$Line <- factor(x = plotdata_sum_df$Line, levels = lines_plot)
# 
# p <- ggplot()
# p <- p + geom_col(data = plotdata_sum_df, mapping = aes(x = Line, y = y_plot), position=position_dodge())
# p <- p + geom_errorbar(data = plotdata_sum_df, mapping = aes(x = Line, y = y_plot, ymin=y_plot-sd_plot, ymax=y_plot+sd_plot), width=.2, position=position_dodge(.9))
# p

