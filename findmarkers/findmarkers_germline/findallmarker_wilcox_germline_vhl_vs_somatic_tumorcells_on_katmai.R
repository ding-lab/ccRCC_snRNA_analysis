# Yige Wu @WashU Apr 2020
## finding differentially expressed gene for each cell type using integrared object

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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/integration/30_aliquot_integration/docker_run_integration/20200212.v3/30_aliquot_integration.20200212.v3.RDS"
srat <- readRDS(file = path_rds)
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)
## set cell type to compare

# make column split -------------------------------------------------------
idaliquot_vhl_germline <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL.Germline) & (bulk_sn_omicsprofile_df$Mut.VHL.Germline != "None")]
idaliquot_vhl_germline
idaliquot_3pdel_only <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL) & (bulk_sn_omicsprofile_df$Mut.VHL == "None")]
idaliquot_vhl_somatic <- as.vector(bulk_sn_omicsprofile_df$Aliquot.snRNA)
idaliquot_vhl_somatic <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL) & (bulk_sn_omicsprofile_df$Mut.VHL != "None") & (bulk_sn_omicsprofile_df$Mut.VHL.Germline == "None")]

# set ident ---------------------------------------------------------------
## make new meta data with marker analysis groups annotated
metadata_new_df <- barcode2celltype_df %>%
  mutate(group_findmarkers = ifelse(Cell_type.shorter == "Tumor cells", 
                                    ifelse(orig.ident %in% idaliquot_vhl_germline, "group1", 
                                           ifelse(orig.ident %in% idaliquot_vhl_somatic, 
                                                  "group2",
                                                  "notincluded")),
                                    "notincluded"))
srat@meta.data <- metadata_new_df
rownames(srat@meta.data) <- metadata_new_df$integrated_barcode
Idents(srat) <- "group_findmarkers"

# run findallmarkers ------------------------------------------------------
markers_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "group1", ident.2 = "group2")
markers_df$row_name <- rownames(markers_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Germline_vs_Somatic_VHL.", "TumorCells.", "FindMarkers.", "Wilcox.", run_id, ".tsv")
write.table(x = markers_df, file = file2write, sep = "\t", quote = F, row.names = F)


