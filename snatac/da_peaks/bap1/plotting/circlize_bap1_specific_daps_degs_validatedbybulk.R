## reference: https://www.cell.com/cell-stem-cell/pdf/S1934-5909(19)30166-3.pdf
## Fig. 6b

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(circlize)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
daps_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_bap1_specific_daps/20210615.v1/BAP1_DAP2Gene.EnhancerPromoter.20210615.v1.tsv")
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/annotate_degs/annotate_bap1_degs_deps_to_chr_regions/20210616.v1/BAP1_vs_PBRM1_NonMutants.DEGs.Chromosome_Regions_Annotated.20210616.v1.tsv")
## input dap deg overlap
dap2deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_bap1_specific_enhancer_promoter_peaks_with_degs/20210615.v1/BAP1_DAP2DEG.20210615.v1.tsv")

# prepare plot data for DAPs-------------------------------------------------------
daps_bed1 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Down") %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC < -1.5, -1.5, avg_log2FC)) %>%
  select(chr, start, end, value)

daps_bed2 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Up") %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC > 1.5, 1.5, avg_log2FC)) %>%
  select(chr, start, end, value)
daps_bed_list = list(daps_bed1, daps_bed2)

# prepare plot data for snRNA DEGs-------------------------------------------------------
degs_snrna_bed1 <- degs_df %>%
  filter(!is.na(Num_sig_up) & (Num_sig_down >=5 & Num_up == 0)) %>%
  filter(!is.na(FDR.cnvcorrected) & FDR.cnvcorrected < 0.05) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_lo2FC.alltumorcells) %>%
  select(chr, start, end, value)
degs_snrna_bed2 <- degs_df %>%
  filter(!is.na(Num_sig_up) & (Num_sig_up >=5 & Num_down == 0)) %>%
  filter(!is.na(FDR.cnvcorrected) & FDR.cnvcorrected < 0.05) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_lo2FC.alltumorcells) %>%
  select(chr, start, end, value)
degs_snrna_bed_list = list(degs_snrna_bed1, degs_snrna_bed2)

# prepare bulk data validating DAP & DEG overlap --------------------------
dap2deg_overlap_df <- dap2deg_df %>%
  filter(DAP_direction == BAP1_vs_OtherTumor_snRNA)

# prepare plot data for bulk RNA-------------------------------------------------------
# > summary(degs_bulkrna_bed2$value)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2079  0.5575  0.9642  1.1671  1.4085  6.0234 
# > summary(degs_bulkrna_bed1$value)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -9.1834 -2.3345 -1.4586 -1.8800 -0.8651 -0.3906 
degs_bulkrna_bed1 <- degs_df %>%
  filter(genesymbol_deg %in% dap2deg_overlap_df$Gene) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < 0.05) %>%
  dplyr::filter(logFC.bulkRNA < 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA < -5, -5, logFC.bulkRNA)) %>%
  select(chr, start, end, value)
degs_bulkrna_bed2 <- degs_df %>%
  filter(genesymbol_deg %in% dap2deg_overlap_df$Gene) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < 0.05) %>%
  filter(logFC.bulkRNA > 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA > 5, 5, logFC.bulkRNA)) %>%
  select(chr, start, end, value)
degs_bulkrna_bed_list = list(degs_bulkrna_bed1, degs_bulkrna_bed2)

# prepare plot data for bulk protein-------------------------------------------------------
deps_bed1 <- degs_df %>%
  filter(genesymbol_deg %in% dap2deg_overlap_df$Gene) %>%
  filter(FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro < 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro < -0.5, -0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value)
deps_bed2 <- degs_df %>%
  filter(genesymbol_deg %in% dap2deg_overlap_df$Gene) %>%
  filter(FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro > 0) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro > 0.5, 0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value)
deps_bed_list = list(deps_bed1, deps_bed2)



# plot version 1 --------------------------------------------------------------------
## initialize
file2write <- paste0(dir_out, "BAP1_specific_daps_degs.bulkRNAFDR0.05.pdf")
pdf(file2write, width = 6, height = 6, useDingbats = F)
circos.par(gap.degree = c(rep(1, 23), 40), start.degree = 90)
circos.initializeWithIdeogram(plotType = c("ideogram", "labels"), labels.cex = 1.2)
# circos.initializeWithIdeogram(plotType = c("axis", "labels"))
## plot track
circos.genomicTrack(daps_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.2, bg.col = "#8DD3C7", bg.border = "grey70")
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
circos.genomicTrack(degs_snrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.2, bg.col = "#FDB462", bg.border = "grey70")
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
circos.genomicTrack(degs_bulkrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.15, bg.col = "#FFFFB3", bg.border = "grey70")
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-5, 0, 5))
circos.genomicTrack(deps_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.15, bg.col = "#BEBADA", bg.border = "grey70")
# circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-0.5, 0, 0.5))
circos.clear()
dev.off()


# plot legend -------------------------------------------------------------
library(ComplexHeatmap)

lgd_points = Legend(at = c("down-regulated in BAP1-mutants", "up-regulated in BAP1-mutants"), type = "points", 
                    legend_gp = gpar(col = c("#377EB8", "#E41A1C"), fill = NA),
                    title = "Point colors")
lgd_background = Legend(at = c("snATAC-seq", "snRNA-seq", "bulk RNA-seq", "bulk proteomics"), type = "grid", 
                    legend_gp = gpar(fill = c("#8DD3C7", "#FDB462", "#FFFFB3", "#BEBADA")),
                    title = "Data types", nrow = 2, border = "grey70")
lgd_list_vertical <- packLegend(lgd_points, lgd_background)
file2write <- paste0(dir_out, "legend.pdf")
pdf(file2write, width = 3, height = 1.5, useDingbats = F)
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
dev.off()
