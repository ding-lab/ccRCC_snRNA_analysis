## reference: https://www.cell.com/cell-stem-cell/pdf/S1934-5909(19)30166-3.pdf

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
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
daps_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_peaks/annotate_BAP1_vs_NonMutant_daps_28samples/20211011.v1/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20211011.v1.tsv")
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/annotate_degs/annotate_bap1_vs_nonmutant_degs_deps_to_chr_regions/20211011.v1/BAP1_vs_NonMutants.DEGs.Chromosome_Regions_Annotated.20211011.v1.tsv")
## input dap deg overlap
dap2deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_degs/overlap_bap1_vs_nonmutant_enhancer_promoter_peaks_28samples_with_degs/20211011.v1/BAP1_vs_NonMutant_DAP2DEG.20211011.v1.tsv")

# set parameters ----------------------------------------------------------
cutoff_log2FC_snATAC <- 0.3
cutoff_log2FC_snRNA <- 0.3
cutoff_fdr_bulkRNA <- 0.0001
cutoff_log2FC_bulkRNA <- 1

# prepare plot data for DAPs-------------------------------------------------------
daps_bed1 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Down") %>%
  # filter(avg_log2FC <= -cutoff_log2FC_snATAC) %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC < -1.5, -1.5, avg_log2FC)) %>%
  select(chr, start, end, value)

daps_bed2 <- daps_df %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction == "Up") %>%
  # filter(avg_log2FC >= cutoff_log2FC_snATAC) %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1]) %>%
  mutate(start = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,2])) %>%
  mutate(end = as.numeric(str_split_fixed(string = peak, pattern = "\\-", n = 3)[,3])) %>%
  mutate(value = ifelse(avg_log2FC > 1.5, 1.5, avg_log2FC)) %>%
  select(chr, start, end, value)
daps_bed_list = list(daps_bed1, daps_bed2)

# prepare plot data for snRNA DEGs-------------------------------------------------------
degs_snrna_bed1 <- degs_df %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  # filter(avg_log2FC.snRNA < -cutoff_log2FC_snRNA) %>%
  filter(Num_sig_down.snRNA >=5) %>%
  # filter(Num_up == 0) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_log2FC.snRNA) %>%
  select(chr, start, end, value)
degs_snrna_bed2 <- degs_df %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  # filter(avg_log2FC.snRNA > cutoff_log2FC_snRNA) %>%
  filter(Num_sig_up.snRNA >=5) %>%
  # filter(Num_up == 0) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  rename(value = avg_log2FC.snRNA) %>%
  select(chr, start, end, value)
degs_snrna_bed_list = list(degs_snrna_bed1, degs_snrna_bed2)

# prepare plot data for bulk protein-------------------------------------------------------
deps_bed1 <- degs_df %>%
  # filter(!is.na(Num_sig_up.snRNA) & (Num_sig_down.snRNA >=5 & Num_up == 0)) %>%
  # filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro < 0) %>%
  # filter(meddiff_exp.bulkpro < -0.1) %>%
  # filter(meddiff_exp.bulkpro < -log2(1.1)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro < -0.5, -0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value, genesymbol_deg)
deps_bed2 <- degs_df %>%
  # filter(!is.na(Num_sig_up.snRNA) & (Num_sig_up.snRNA >=5 & Num_down == 0)) %>%
  # filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(FDR.bulkpro < 0.05) %>%
  filter(meddiff_exp.bulkpro > 0) %>%
  # filter(meddiff_exp.bulkpro > 0.1) %>%
  # filter(meddiff_exp.bulkpro > log2(1.1)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(meddiff_exp.bulkpro > 0.5, 0.5, meddiff_exp.bulkpro)) %>%
  select(chr, start, end, value, genesymbol_deg)
deps_bed_list = list(deps_bed1, deps_bed2)

# prepare plot data for bulk RNA-------------------------------------------------------
degs_bulkrna_bed1 <- degs_df %>%
  # filter(!is.na(Num_sig_up.snRNA) & (Num_sig_down.snRNA >=5 & Num_up == 0)) %>%
  # filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter(logFC.bulkRNA < -(cutoff_log2FC_bulkRNA)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA < -5, -5, logFC.bulkRNA)) %>%
  select(chr, start, end, value, genesymbol_deg)
degs_bulkrna_bed2 <- degs_df %>%
  # filter(!is.na(Num_sig_up.snRNA) & (Num_sig_up.snRNA >=5 & Num_down == 0)) %>%
  # filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter(logFC.bulkRNA > cutoff_log2FC_bulkRNA) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = ifelse(logFC.bulkRNA > 5, 5, logFC.bulkRNA)) %>%
  select(chr, start, end, value, genesymbol_deg)
degs_bulkrna_bed_list = list(degs_bulkrna_bed1, degs_bulkrna_bed2)

# prepare plot data for highlighting genes-------------------------------------------------------
# genes_highlight <- dap2deg_df %>%
#   filter(abs(avg_log2FC.snATAC) >= 0.5 & abs(avg_log2FC.snRNA) >= 0.5)
genes_highlight <- dap2deg_df %>%
  filter(abs(avg_log2FC.snATAC) >= cutoff_log2FC_snATAC & abs(avg_log2FC.snRNA) >= cutoff_log2FC_snRNA) %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA) & (avg_log2FC.snATAC >= 0 & avg_log2FC.snRNA >= 0 | avg_log2FC.snATAC <= 0 & avg_log2FC.snRNA <= 0))
genes_validatedbybulkRNA_df <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkRNA) & FDR.bulkRNA < cutoff_fdr_bulkRNA) %>%
  filter((Num_sig_up.snRNA >=5 & logFC.bulkRNA > 0) | (Num_sig_down.snRNA >=5 & logFC.bulkRNA < 0)) %>%
  filter(!is.na(start_position)) %>%
  mutate(chr = paste0("chr", chromosome_name))
genes_validatedbybulkpro_df <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(FDR.snRNA.cnvcorrected) & FDR.snRNA.cnvcorrected < 0.05) %>%
  filter(!is.na(FDR.bulkpro) & FDR.bulkpro < 0.05) %>%
  filter((Num_sig_up.snRNA >=5 & meddiff_exp.bulkpro > 0) | (Num_sig_down.snRNA >=5 & meddiff_exp.bulkpro < 0)) %>%
  filter(!is.na(start_position))

highlight_bed <- degs_df %>%
  filter(genesymbol_deg %in% genes_highlight$Gene) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(start = start_position) %>%
  rename(end = end_position) %>%
  mutate(value = genesymbol_deg) %>%
  mutate(validated_by_bulkRNA = (genesymbol_deg %in% genes_validatedbybulkRNA_df$genesymbol_deg)) %>%
  mutate(validated_by_bulkprotein = (genesymbol_deg %in% genes_validatedbybulkpro_df$genesymbol_deg)) %>%
  mutate(deg_group = ifelse(validated_by_bulkRNA, 
                            ifelse(validated_by_bulkprotein, "bulk RNA & bulk protein", "bulk RNA only "),
                            ifelse(validated_by_bulkprotein, "bulk protein only", "not validated by bulk"))) %>%
  # mutate(deg_color = ifelse(validated_by_bulkRNA, 
  #                           ifelse(validated_by_bulkprotein, "#984EA3", "#D95F02"),
  #                           ifelse(validated_by_bulkprotein, "#E7298A", "black"))) %>%
  mutate(deg_color = ifelse(validated_by_bulkRNA | validated_by_bulkprotein, 
                            ifelse(avg_log2FC.snRNA > 0, "#E41A1C" , "#377EB8"),
                            ifelse(avg_log2FC.snRNA > 0, "#FB9A99" , "#A6CEE3"))) %>%
  mutate(deg_font = ifelse(!validated_by_bulkRNA & !validated_by_bulkprotein, 3, 4))  %>%
  select(chr, start, end, value, deg_group, deg_color, deg_font)

# plot all data tracks --------------------------------------------------------------------
colors_datatypes <- c("#99D8C9", "#FFFFB3", "#FDD0A2",  "#DECBE4")
names(colors_datatypes) <- c("snATAC", "snRNA", "bulkRNA", "bulkprotein")
## initialize
file2write <- paste0(dir_out, "BAP1_specific_daps_degs.bulkRNAFDR", cutoff_fdr_bulkRNA, 
                     ".log2FC",  cutoff_log2FC_bulkRNA,
                     ".log2snATAC", cutoff_log2FC_snATAC,
                     ".log2snRNA", cutoff_log2FC_snRNA,
                     ".pdf")
pdf(file2write, width = 7, height = 7, useDingbats = F)
# pdf(file2write, width = 7.5, height = 7.5, useDingbats = F)
circos.par(gap.degree = c(rep(1, 23), 40), start.degree = 90)
circos.initializeWithIdeogram(plotType = NULL)
# circos.initializeWithIdeogram(plotType = c("labels"), track.height = 0.05)
## plot labels
### font: 1=regular, 2=bold, 3=italic, 4=bold-italic
circos.genomicLabels(highlight_bed, labels.column = 4, side = "outside", connection_height = 0.08, labels_height = 0.08, cex = 0.8, col = highlight_bed$deg_color, line_lwd = 0.5, font = highlight_bed$deg_font)
# circos.genomicIdeogram()
## plot track
circos.genomicTrack(daps_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["snATAC"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
circos.genomicTrack(degs_snrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["snRNA"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
# circos.genomicLabels(highlight_bed, labels.column = 4, side = "inside", connection_height = 0.01, labels_height = 0.03, cex = 0.5)
circos.genomicTrack(degs_bulkrna_bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["bulkRNA"], bg.border = NA)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-5, 0, 5))

circos.genomicTrack(deps_bed_list,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      color_i <- c("#377EB8", "#E41A1C")[i]
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
                    }, track.height = 0.14, bg.col = colors_datatypes["bulkprotein"], bg.border = NA)
# circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8)
circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-0.5, 0, 0.5))
circos.clear()
dev.off()

# plot version 1 --------------------------------------------------------------------
## initialize
# file2write <- paste0(dir_out, "BAP1_specific_daps_degs.onlysndata.pdf")
# pdf(file2write, width = 6, height = 6, useDingbats = F)
# circos.par(gap.degree = c(rep(1, 23), 40), start.degree = 90)
# circos.initializeWithIdeogram(plotType = NULL)
# # circos.initializeWithIdeogram(plotType = c("axis", "labels"))
# ## plot labels
# circos.genomicLabels(highlight_bed, labels.column = 4, side = "outside", connection_height = 0.01, labels_height = 0.03, cex = 0.5, col = highlight_bed$deg_group)
# circos.genomicIdeogram()
# ## plot track
# circos.genomicTrack(daps_bed_list, 
#                     panel.fun = function(region, value, ...) {
#                       i = getI(...)
#                       color_i <- c("#377EB8", "#E41A1C")[i]
#                       circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
#                     }, track.height = 0.15, bg.col = "#8DD3C7", bg.border = "grey70")
# circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
# circos.genomicTrack(degs_snrna_bed_list, 
#                     panel.fun = function(region, value, ...) {
#                       i = getI(...)
#                       color_i <- c("#377EB8", "#E41A1C")[i]
#                       circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = color_i, ...)
#                     }, track.height = 0.15, bg.col = "#FDB462", bg.border = "grey70")
# circos.yaxis(side = "left", sector.index = "chr1", labels.cex = 0.8, at = c(-1.5, 0, 1.5))
# circos.clear()
# dev.off()

