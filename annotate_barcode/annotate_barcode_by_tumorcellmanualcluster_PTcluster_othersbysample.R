# Yige Wu @WashU Aug 2021

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
## input cell type annotation
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv")
## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
## input PT + LOH reclustered
pt_loh_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_PT_LOH_merged_katmai/20210809.v1/PT_LOH_Merged.Metadata.20210809.v1.tsv")
## input id data
id_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")


# process PT & LOH barcodes -----------------------------------------------------
## process PT barcodes
ptloh_bycluster_df <- pt_loh_df %>%
  mutate(cell_group = paste0(ifelse(seurat_clusters %in% c(9, 4, 11, 3, 0, 8, 15, 13, 5), "PT_C", "LOH_C"), seurat_clusters))
barcodes_pt_loh <- ptloh_bycluster_df$barcode

# annotate tumor barcodes -------------------------------------------------------
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(sample_barcode = paste0(orig.ident, "_", barcode))
tumor_bycluster_df <- barcode2celltype_df %>%
  mutate(sample_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  filter(sample_barcode %in% barcode2tumorsubcluster_df$sample_barcode)
tumor_bycluster_df$cell_group <- mapvalues(x = tumor_bycluster_df$sample_barcode, from = barcode2tumorsubcluster_df$sample_barcode, to = as.vector(barcode2tumorsubcluster_df$Cluster_Name))
barcodes_tumor <- tumor_bycluster_df$sample_barcode

# process other barcodes --------------------------------------------------
barcodes_other_df <- barcode2celltype_df %>%
  mutate(sample_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(cell_group = "other") %>%
  select(sample_barcode, cell_group)

# combine -----------------------------------------------------------------
cellgroup_df <- rbind(ptloh_bycluster_df %>%
                        mutate(sample_barcode = barcode) %>%
                        select(sample_barcode, cell_group),
                      tumor_bycluster_df %>%
                        select(sample_barcode, cell_group),
                      barcodes_other_df %>%
                        select(sample_barcode, cell_group))
table(cellgroup_df$cell_group)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode_TumorPTLOH_Reclustered.", run_id, ".tsv")
write.table(x = cellgroup_df, file = file2write, sep = "\t", row.names = F, quote = F)
