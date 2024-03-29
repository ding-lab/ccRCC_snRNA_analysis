# Yige Wu @WashU Sep 2021
## source activate ccrcc_snrna
## run ulimit -s unlimited

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
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
## library additional libaries
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
library(Signac)
library(ggplot2)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 7 * 1024^3) # for 7 Gb RAM

# input dependencies ------------------------------------------------------
## input the merged object
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v3.20210725.rds',sep=''))
## input barcode to tumor clusters
# barcode2cluster_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/RNA_clusters_Tumor_PT_clusters.snATAC.20210917.tsv")
barcode2cluster_df <- fread(data.table = F, input = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/12.Transfer_tumorClusters_fromRNA/TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters_Tumor_PT_clusters.snATAC.20210917.tsv")
## input tumor cluster grouping
cluster2group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_epithelial_group_bytumorcluster/20210927.v1/Tumorcluster_EpithelialGroup.20210927.v1.tsv")
## input the peaks to plot
peak2deg_motif_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/filter_motifs/filter_motifs_for_EMT_vs_selectedEpithelial_peaks_for_degs/20210928.v1/Epithelial_up_DEG_DAP_DAM_Overlap.20210928.v1.tsv")

# preprocess ATAC object --------------------------------------------------
print("Start mapping cell_group")
atac@meta.data$cell_group <- mapvalues(x = rownames(atac@meta.data), from = barcode2cluster_df$sample_barcode, to = as.vector(barcode2cluster_df$cell_group), warn_missing = F)
atac@meta.data$cell_group[atac@meta.data$cell_group == rownames(atac@meta.data)] <- "other"
table(atac@meta.data$cell_group)
print("Finished mapping cell_group")
print("Start mapping epithelial group")
cluster2group_df <- cluster2group_df %>%
  mutate(cluster_name.formatted = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
cluster2group_df <- merge(x = data.frame(cluster_name.formatted = unique(barcode2cluster_df$cell_group)), y = cluster2group_df, by = c("cluster_name.formatted"), all.x = T)
cluster2group_df$epithelial_group[is.na(cluster2group_df$epithelial_group)] <- "other"
atac@meta.data$epithelial_group <- mapvalues(x = atac@meta.data$cell_group, from = cluster2group_df$cluster_name.formatted, to = as.vector(cluster2group_df$epithelial_group), warn_missing = F)
atac@meta.data$epithelial_group[atac@meta.data$cell_group == atac@meta.data$epithelial_group] <- "other"
table(atac@meta.data$epithelial_group)
print("Finished mapping epithelial_group")
cluster2group_df <- cluster2group_df %>%
  arrange(factor(epithelial_group, levels = c("EMT", "Epithelial-weak", "Epithellal-intermediate", "Epithelial-strong", "other")), score)
atac@meta.data$cell_group=factor(atac@meta.data$cell_group, levels=c(cluster2group_df$cluster_name.formatted, "other"))
print("Start changing ident")
Idents(atac) <- "cell_group"
print("Finished changing ident")

print("Start subsetting")
atac_subset=subset(atac,(cell_group %in% c("C3L-00079-T1_C4", "C3L-01302-T1_C1",
                                           "C3L-00088-T2_C1", "C3N-00733-T1_C1", "C3L-00416-T2_C1", "C3L-00010-T1_C1", "C3L-00088-T1_C1")))
# atac_subset=subset(atac,!(cell_group %in% c("other", cluster2group_df$cluster_name.formatted[cluster2group_df$epithelial_group == "other"])))
# atac_subset=atac
print("Finished subsetting")
rm(atac)

# make colors -------------------------------------------------------------
## make colors
colors_tumorgroup_sim <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(1, 5, 4, 2, 9)]
names(colors_tumorgroup_sim) <- c("EMT",  "Epithelial-weak", "Epithellal-intermediate", "Epithelial-strong", "other")
colors_tumorgroup <- colors_tumorgroup_sim[cluster2group_df$epithelial_group]
names(colors_tumorgroup) <- cluster2group_df$cluster_name.formatted

# process peaks -----------------------------------------------------------
plotdata_df <- peak2deg_motif_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  filter(Gene == "EPB41L4A") %>%
  filter(avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1) %>%
  filter(!is.na(avg_log2FC.tf.snRNA) & avg_log2FC.tf.snRNA < 0) %>%
  # filter(motif.name %in% c("JUN", "TWIST1")) %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak, motif.name, motif_coord, avg_log2FC.tf.snRNA) %>%
  unique()
print("Finished processing peaks")


# plot --------------------------------------------------------------------
peak_plot <- unique(plotdata_df$peak)
print(paste0("Start processing peak: ", peak_plot))

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
gene_plot <- plotdata_df$Gene[plotdata_df$peak == peak_plot]
motif_coord1 <- plotdata_df$motif_coord[plotdata_df$motif.name == "NFIB"]
motif_coord2 <- plotdata_df$motif_coord[plotdata_df$motif.name == "NFIA"]
motif_coords <- plotdata_df$motif_coord[plotdata_df$motif.name %in% c("NFIB", "NFIA")]

# plot --------------------------------------------------------------------
print(paste0("Start processing cov_plot"))
cov_obj= Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded,
  annotation = F, 
  peaks = F,
  links=FALSE)
cov_obj <- cov_obj + scale_fill_manual(values =  colors_tumorgroup)
# cov_obj <- cov_obj +  theme_classic(base_size = 12) + theme(legend.position = "none")
print("Finished cov_plot")

peak_plot_obj <- Signac::PeakPlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  peaks = StringToGRanges(peak_plot, sep = c("-", "-")))
# peak_plot_obj <- peak_plot_obj  + theme_classic(base_size = 12)
print("Finished peak plot")

gene_plot_obj <- Signac::AnnotationPlot(
  object = atac_subset,
  region = peak_plot_expanded)
# gene_plot_obj <- gene_plot_obj  + theme_classic(base_size = 12)

# motif_plot_obj1 <- Signac::PeakPlot(
#   object = atac_subset,
#   region = peak_plot_expanded, 
#   peaks = StringToGRanges(motif_coord1, sep = c("-", "-")))
# 
# motif_plot_obj2 <- Signac::PeakPlot(
#   object = atac_subset,
#   region = peak_plot_expanded, 
#   peaks = StringToGRanges(motif_coord2, sep = c("-", "-")))

motif_plot_obj <- Signac::PeakPlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  peaks = StringToGRanges(motif_coords, sep = c("-", "-")))

# p <- Signac::CombineTracks(
#   plotlist = list(cov_obj, peak_plot_obj, motif_plot_obj1, motif_plot_obj2, gene_plot_obj),
#   heights = c(6, 0.3, 0.3, 0.3, 1))

p <- Signac::CombineTracks(
  plotlist = list(cov_obj, peak_plot_obj,motif_plot_obj, gene_plot_obj),
  heights = c(6, 0.3, 1))

print("Finished CombineTracks")

## write output
# file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, gene_plot, "_", gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".pdf")
pdf(file2write, width = 6, height = 6)
print(p)
dev.off()

