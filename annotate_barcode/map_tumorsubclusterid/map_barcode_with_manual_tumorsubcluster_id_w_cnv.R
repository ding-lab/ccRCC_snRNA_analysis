# Yige Wu @WashU Jun 2020
## annotate each barcode to the manual tumor subcluster number (raw)
## annotate cells with mixture marker gene expression to be unknown
## annotate cells with hallmark cnv abeit with mixture marker gene expression to be tumor cells

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
## input barcode2seurat cluster info
barcode2seuratcluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_seurat_tumorsubcluster_id/20200603.v1/Barcode2SeuratClusterID.20200603.v1.tsv")
## input the cell to cell type table
sratcluster2manualcluster_df <- readxl::read_xlsx(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Tumor_Subcluster/Individual.TumorCluster2Cell_Type.20200616.v1.xlsx")
## input barcode to cluster mapping table from 
all_integrated_barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input cnv info by cell for rescueing unknown cells
cnv_state_bycell_bygene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/annotate_barcode_with_gene_level_cnv_using_cnv_genes/20200518.v1/Individual.20200305.v1.CNV_State_By_Gene_By_Barcode.20200518.v1.tsv")
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# get barcode2manualcluster -----------------------------------------------
barcode2manualcluster_df <- merge(x = barcode2seuratcluster_df %>%
                                    select(orig.ident, id_seurat_cluster, barcode), 
                                  y = sratcluster2manualcluster_df, 
                                  by.x = c("orig.ident", "id_seurat_cluster"), 
                                  by.y = c("Aliquot", "id_seurat_cluster"), all.x = T)

# divide barcodes to those from certain tumor cell cluster and those with mixture marker gene expression ---------------------------
barcode2manualcluster_df1 <- barcode2manualcluster_df %>%
  filter(Is_Malignant != "Mixture") %>%
  rename(individual_barcode = barcode) %>%
  rename(Id_TumorManualCluster = id_manual_cluster) %>%
  mutate(Cell_type.shorter = "Tumor cells") %>%
  mutate(Cell_type.detailed = Cell_type.shorter) %>%
  rename(Most_Enriched_Cell_Group = Enriched_Cell_Group) %>%
  rename(Most_Enriched_Cell_Type1 = Enriched_Cell_Type1) %>%
  rename(Most_Enriched_Cell_Type2 = Enriched_Cell_Type2) %>%
  rename(Most_Enriched_Cell_Type3 = Enriched_Cell_Type3) %>%
  rename(Most_Enriched_Cell_Type4 = Enriched_Cell_Type4) %>%
  mutate(Most_Enriched_Cell_Type2 = "") %>%
  mutate(Most_Enriched_Cell_Type3 = "") %>%
  mutate(Most_Enriched_Cell_Type4 = "") %>%
  select(orig.ident, individual_barcode,
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)

# get barcodes with hallmark cnvs and change the cell type of those within the mixture marker gene expression-----------------------------------------
genes_hallmarkcnvs <- knowncnvgenes_df$Gene_Symbol[grepl(x = knowncnvgenes_df$Cytoband, pattern = "3p|5q|14q")]
rescue_cnv_state_bycell_bygene_df <- cnv_state_bycell_bygene_df %>%
  filter(gene_expected_cna_state == cna_state) %>%
  filter(gene_symbol %in% genes_hallmarkcnvs) %>%
  select(id_aliquot, barcode_individual) %>%
  unique() %>%
  mutate(has_hallmarkcnv = T)


# label barcodes with mixture marker gene expression ----------------------
barcode2manualcluster_df2 <- barcode2manualcluster_df %>%
  filter(Is_Malignant == "Mixture")
barcode2manualcluster_df2 <- merge(barcode2manualcluster_df2, rescue_cnv_state_bycell_bygene_df,
                                   by.x = c("orig.ident", "barcode"), 
                                   by.y = c("id_aliquot", "barcode_individual"), all.x = T)
table(barcode2manualcluster_df2$has_hallmarkcnv)
which(is.na(barcode2manualcluster_df2$has_hallmarkcnv))
## label tumor cells
barcode2manualcluster_df2$Enriched_Cell_Group[!is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- "Nephron_Epithelium"
barcode2manualcluster_df2$Enriched_Cell_Type1[!is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <-"Proximal tubule"
barcode2manualcluster_df2$Enriched_Cell_Type2[!is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
barcode2manualcluster_df2$Enriched_Cell_Type3[!is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
barcode2manualcluster_df2$Enriched_Cell_Type4[!is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
## label unknown
barcode2manualcluster_df2$Enriched_Cell_Group[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- "Unknown"
barcode2manualcluster_df2$Enriched_Cell_Type1[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <-""
barcode2manualcluster_df2$Enriched_Cell_Type2[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
barcode2manualcluster_df2$Enriched_Cell_Type3[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
barcode2manualcluster_df2$Enriched_Cell_Type4[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- ""
barcode2manualcluster_df2$id_manual_cluster[is.na(barcode2manualcluster_df2$has_hallmarkcnv)] <- NA
## format
barcode2manualcluster_df2 <- barcode2manualcluster_df2 %>%
  rename(individual_barcode = barcode) %>%
  rename(Id_TumorManualCluster = id_manual_cluster) %>%
  mutate(Cell_type.shorter = ifelse(Enriched_Cell_Group == "Nephron_Epithelium", "Tumor cells", "Unknown")) %>%
  mutate(Cell_type.detailed = Cell_type.shorter) %>%
  rename(Most_Enriched_Cell_Group = Enriched_Cell_Group) %>%
  rename(Most_Enriched_Cell_Type1 = Enriched_Cell_Type1) %>%
  rename(Most_Enriched_Cell_Type2 = Enriched_Cell_Type2) %>%
  rename(Most_Enriched_Cell_Type3 = Enriched_Cell_Type3) %>%
  rename(Most_Enriched_Cell_Type4 = Enriched_Cell_Type4) %>%
  select(orig.ident, individual_barcode,
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)

# merge -------------------------------------------------------------------
tumor_barcode2celltype_df <- rbind(barcode2manualcluster_df1, barcode2manualcluster_df2)
nrow(tumor_barcode2celltype_df)

## add integrated barcode
tumor_barcode2celltype_df <- merge(tumor_barcode2celltype_df, 
                                  all_integrated_barcode2cluster_df %>%
                                    rename(integrated_barcode = barcode) %>%
                                    mutate(individual_barcode = str_split_fixed(string = integrated_barcode, pattern = "_", n = 2)[,1]) %>%
                                    select(orig.ident, individual_barcode, integrated_barcode),
                                  by = c("orig.ident", "individual_barcode"),
                                  all.x = T)
nrow(tumor_barcode2celltype_df)
## [1] 93277
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2TumorSubclusterId.", run_id, ".tsv")
write.table(x = tumor_barcode2celltype_df, file = file2write, sep = '\t', quote = F, row.names = F)


