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
## input PT + LOH reclustered
pt_loh_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_PT_LOH_merged_katmai/20210809.v1/PT_LOH_Merged.Metadata.20210809.v1.tsv")
## input id data
id_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# annotate barcodes -------------------------------------------------------
table(barcode2celltype_df$Cell_group_w_epithelialcelltypes)
table(barcode2celltype_df$Cell_group14_w_transitional)
## process PT barcodes
pt_bycluster_df <- pt_loh_df %>%
  filter(seurat_clusters %in% c(9, 4, 11, 3, 0, 8, 15, 13, 5)) %>%
  mutate(cell_group = paste0("PT_C", seurat_clusters))
barcodes_pt <- pt_bycluster_df$barcode
## process non-pt barcodes
barcode2celltype_df$easyid <- mapvalues(x = barcode2celltype_df$orig.ident, from = id_df$Aliquot.snRNA, to = as.vector(id_df$Aliquot.snRNA.WU))
nonpt_df <- barcode2celltype_df %>%
  mutate(barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  filter(!(barcode %in% barcodes_pt)) %>%
  mutate(cell_group = paste0(ifelse(Cell_group14_w_transitional == "EMT tumor cells", "EMT tumor cells", Cell_group_w_epithelialcelltypes), "_", easyid))
## combine
cellgroup_df <- rbind(pt_bycluster_df %>%
                        select(barcode, cell_group),
                      nonpt_df %>%
                        select(barcode, cell_group))
table(cellgroup_df$cell_group)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode_byEpithelialCelltypes_BySample.", run_id, ".tsv")
write.table(x = cellgroup_df, file = file2write, sep = "\t", row.names = F, quote = F)
