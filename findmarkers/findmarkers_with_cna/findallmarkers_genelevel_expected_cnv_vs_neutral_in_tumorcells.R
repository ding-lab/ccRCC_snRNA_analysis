# Yige Wu @WashU May 2020
## run DEG analysis for cells with certain CNV vs CN neutral

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for individual seurat object
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/annotate_barcode_with_gene_level_cnv_using_cnv_genes/20200518.v1/Individual.20200305.v1.CNV_State_By_Gene_By_Barcode.20200518.v1.tsv", data.table = F)

# for each aliquot, input seurat object and fetch data and write data --------------------
## initiate the vector to record the file path
aliquot_vec <- NULL
gene_symbol_vec <- NULL
deg_filepath_vec <- NULL
markers_wilcox_df <- NULL
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
  ## make output directory for this aliquot
  dir_out_sub1 <- paste0(dir_out, aliquot_tmp, "/")
  dir.create(dir_out_sub1)
  
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
    
  ## identify genes with detected cna state
  genes2process <- barcode2cnv_df %>%
    filter(id_aliquot == aliquot_tmp) %>%
    select(gene_symbol) %>%
    unique()
  genes2process <- as.vector(genes2process$gene_symbol)
  
  for (gene_symbol_tmp in genes2process) {
    ## get expected cna_state
    gene_expected_cna_state_tmp <- unique(barcode2cnv_df$gene_expected_cna_state[barcode2cnv_df$gene_symbol == gene_symbol_tmp])
    ## make new meta data with cna state for the current gene
    meta.data_cna <- barcode2cnv_df %>%
      filter(id_aliquot == aliquot_tmp) %>%
      filter(gene_symbol == gene_symbol_tmp)
    rownames(meta.data_cna) <- meta.data_cna$barcode_individual
    ## count the number of cells with expected cna state and neutral
    num_cells_expectedcna <- meta.data_cna %>%
      filter(cna_state == gene_expected_cna_state) %>%
      nrow()
    num_cells_neutralcna <- meta.data_cna %>%
      filter(cna_state == "Neutral") %>%
      nrow()
    if (num_cells_expectedcna <= 10 || num_cells_neutralcna <= 10) {
      next()
    }
    ### set the identities to the cna state of the current gene
    srat_tmp <- srat
    srat_tmp@meta.data <- meta.data_cna
    Idents(srat_tmp) <- "cna_state"
    
    ### find all markers using Wilcox testing
    markers_wilcox <- tryCatch(expr = ,FindMarkers(object = srat_tmp, ident.1 = gene_expected_cna_state_tmp, ident.2 = "Neutral", logfc.threshold = 0.1,
                                                  only.pos = FALSE, test.use = "wilcox"),
                               error = function(e) {warning("Marker identification failed.");return(NULL)})
    if (is.null(markers_wilcox)) {
      next()
    }
    markers_wilcox$de_gene_symbol <- rownames(markers_wilcox)
    markers_wilcox <- markers_wilcox %>%
      mutate(cna_gene_symbol = gene_symbol_tmp) %>%
      mutate(id_aliquot = aliquot_tmp) %>%
      mutate(cna_gene_state.1 = gene_expected_cna_state_tmp) %>%
      mutate(cna_gene_state.2 = "Neutral") %>%
      mutate(num_cells.1 = num_cells_expectedcna) %>%
      mutate(num_cells.2 = num_cells_neutralcna)
    ### write table
    file2write <- paste0(dir_out_sub1, aliquot_tmp, ".", gene_symbol_tmp, ".", "ExpectedCNV_vs_Neutral", ".", ".FindAllMarkers.Wilcox.", ".", run_id, ".tsv")
    write.table(markers_wilcox, file = file2write, quote = F, sep = "\t", row.names = F)
    
    ### store file path
    deg_filepath_vec <- c(file2write, deg_filepath_vec)
    gene_symbol_vec <- c(gene_symbol_tmp, gene_symbol_vec)
    aliquot_vec <- c(aliquot_tmp, aliquot_vec)
    
    ## store deg result
    markers_wilcox_df <- rbind(markers_wilcox_df, markers_wilcox)
  }
}

### write table
file2write <- paste0(dir_out, "ExpectedCNV_vs_Neutral", ".", ".FindAllMarkers.Wilcox.", ".", run_id, ".tsv")
write.table(markers_wilcox_df, file = file2write, quote = F, sep = "\t", row.names = F)
