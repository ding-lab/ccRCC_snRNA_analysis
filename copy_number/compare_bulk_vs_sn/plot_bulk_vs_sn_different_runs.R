# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
### with InferCNV run (20200305.v1)
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200319.v1/bulk_sn_omics_profile.20200319.v1.tsv", data.table = F)
## input InferCNV run (20200207.v1)
frac_expectedcnv_bychr_byaliquot_df <- fread(input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/estimate_fraction_of_tumorcells_with_expectedcnv_perchrregion_per_sample_using_cnvgenes/20200415.v1/Individual.20200207.v1fraction_of_tumorcells.expectedCNA.by_chr_region.20200415.v1.tsv")

# make heatmap body -------------------------------------------------------
plot_data_df <- bulk_sn_omicsprofile_df[, c("CN.bulk.3p", "CN.bulk.5q", "CN.bulk.14q")]
plot_data_df$id_aliquot <- bulk_sn_omicsprofile_df$Aliquot.snRNA
## order the samples by their CNV profile
plot_data_df <- plot_data_df %>%
  mutate(cnv_combined_status = paste0(CN.bulk.3p, "_", CN.bulk.5q, "_", CN.bulk.14q)) %>%
  arrange(cnv_combined_status)
## get the aliquot ids in order
ids_aliquot <- plot_data_df$id_aliquot
ids_aliquot
## make matrix out of data frame
plot_data_mat <- t(as.matrix(plot_data_df[,c("CN.bulk.3p", "CN.bulk.5q", "CN.bulk.14q")]))
plot_data_mat
colnames(plot_data_mat) <- ids_aliquot

# make heatmap body colors ------------------------------------------------
colors_heatmapbody = structure(c("red", "blue", "white"), names = c("gain", "loss", "neutral")) # black, red, green, blue

# make top column annotation --------------------------------------------------
merged_cnv_df <- merge(bulk_sn_omicsprofile_df %>%
                         select(Aliquot.snRNA, CN.sn.3p_loss.fraction, CN.sn.5q_gain.fraction, CN.sn.14q_loss.fraction) %>%
                         rename(CN.sn.3p_loss.fraction.Ref_LabeledSame = CN.sn.3p_loss.fraction) %>%
                         rename(CN.sn.5q_loss.fraction.Ref_LabeledSame = CN.sn.5q_gain.fraction) %>%
                         rename(CN.sn.14q_loss.fraction.Ref_LabeledSame = CN.sn.14q_loss.fraction) %>%
                         rename(aliquot = Aliquot.snRNA),
                       frac_expectedcnv_bychr_byaliquot_df %>%
                         select(aliquot, `3p`, `5q`, `14q`) %>%
                         rename(CN.sn.3p_loss.fraction.Ref_LabeledByCellGroup = `3p`) %>%
                         rename(CN.sn.5q_gain.fraction.Ref_LabeledByCellGroup = `5q`) %>%
                         rename(CN.sn.14q_loss.fraction.Ref_LabeledByCellGroup = `14q`),
                       by = c("aliquot"), all = T)
top_col_anno_df <- merged_cnv_df %>%
  select(-aliquot)
rownames(top_col_anno_df) <- as.vector(merged_cnv_df$aliquot)
top_col_anno_df <- top_col_anno_df[ids_aliquot,]
rownames(top_col_anno_df) <- ids_aliquot
### make neutral variant info for the normal sample
normal_aliquot_ids <- id_metadata_df$Aliquot.snRNA[id_metadata_df$Sample_Type == "Normal"]
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids, paste0("CN.", c('3p', '5q', "14q"))] <- "neutral"
# ### make color for tumor purity
# tumorpurity_color_fun <-  colorRamp2(c(quantile(top_col_anno_df$TumorPurity.snRNA, 0.1, na.rm=T), 
#                                        quantile(top_col_anno_df$TumorPurity.snRNA, 0.5, na.rm=T), 
#                                        quantile(top_col_anno_df$TumorPurity.snRNA, 0.9, na.rm=T)),c("blue", "white", "red"))
# 
## top column annotation object
top_col_anno = HeatmapAnnotation(CN.sn.3p_loss.fraction.Ref_LabeledByCellGroup = anno_barplot(top_col_anno_df$CN.sn.3p_loss.fraction.Ref_LabeledByCellGroup, 
                                                                                              height = unit(5, "mm")),
                                 CN.sn.5q_gain.fraction.Ref_LabeledByCellGroup = anno_barplot(x = top_col_anno_df$CN.sn.5q_gain.fraction.Ref_LabeledByCellGroup, 
                                                                                              height = unit(5, "mm")),
                                 CN.sn.14q_loss.fraction.Ref_LabeledByCellGroup = anno_barplot(x = top_col_anno_df$CN.sn.14q_loss.fraction.Ref_LabeledByCellGroup, 
                                                                                               height = unit(5, "mm")),
                                 CN.sn.3p_loss.fraction.Ref_LabeledSame = anno_barplot(top_col_anno_df$CN.sn.3p_loss.fraction.Ref_LabeledSame, 
                                                                                       height = unit(5, "mm")),
                                 CN.sn.5q_loss.fraction.Ref_LabeledSame = anno_barplot(x = top_col_anno_df$CN.sn.5q_loss.fraction.Ref_LabeledSame, 
                                                                                       height = unit(5, "mm")),
                                 CN.sn.14q_loss.fraction.Ref_LabeledSame = anno_barplot(x = top_col_anno_df$CN.sn.14q_loss.fraction.Ref_LabeledSame, 
                                                                                        height = unit(5, "mm")))



# make the whole heatmap --------------------------------------------------
p <- Heatmap(matrix = plot_data_mat,
             col = colors_heatmapbody,
             na_col = "grey", top_annotation = top_col_anno,
             show_heatmap_legend = T)
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "test", ".png")
png(filename = file2write,width = 1400, height = 500, res = 150)
print(p)
dev.off()
