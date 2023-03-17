## reference: https://www.cell.com/cell-stem-cell/pdf/S1934-5909(19)30166-3.pdf
## Fig. 6b

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
  "ggplot2",
  "circlize"
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
peaks_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_bap1_specific_daps/20210615.v1/BAP1_DAP2Gene.EnhancerPromoter.20210615.v1.tsv")

# set parameters ----------------------------------------------------------
cutoff_log2FC_snATAC <- 0.3
cutoff_log2FC_snRNA <- 0.3

# prepare plot data -------------------------------------------------------
plotdata_bed1 <- peaks_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Down") %>%
  filter(avg_log2FC <= -cutoff_log2FC_snATAC) %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC < -1.5, -1.5, avg_log2FC)) %>%
  select(chr, start, end, value)

plotdata_bed2 <- peaks_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Up") %>%
  filter(avg_log2FC >= cutoff_log2FC_snATAC) %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC > 1.5, 1.5, avg_log2FC)) %>%
  select(chr, start, end, value)

bed_list = list(plotdata_bed1, plotdata_bed2)

# plot version 1 --------------------------------------------------------------------
## initialize
file2write <- paste0(dir_out, "BAP1_specific_peaks.pdf")
pdf(file2write, width = 6, height = 6, useDingbats = F)
circos.par(gap.degree = c(rep(1, 23), 30), start.degree = 90)
circos.initializeWithIdeogram(plotType = c("ideogram", "labels"), labels.cex = 1.2)
# circos.initializeWithIdeogram(plotType = c("axis", "labels"))
## plot track
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.2)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 1.2, at = c(-1.5, 0, 1.5))

## add density
# circos.genomicDensity(bed_list, col = c("#377EB8", "#E41A1C"), track.height = 0.2, overlap = F)
circos.genomicDensity(bed_list[[1]], col = c("#377EB8"), track.height = 0.2)
circos.genomicDensity(bed_list[[2]], col = c("#E41A1C"), track.height = 0.2)
circos.clear()
dev.off()
## save source data
source_data_df <- rbind(plotdata_bed1, plotdata_bed2)
write.table(x = source_data_df, file = "~/Desktop/SF6b.SourceData.tsv", quote = F, sep = "\t", row.names = F)
