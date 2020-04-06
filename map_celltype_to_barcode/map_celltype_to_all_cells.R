# Yige Wu @WashU Apr 2020
## make barcode to cell type mapping table for the integrated dataset

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode to cluster mapping table from 
all_integrated_barcode2cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input cluster to cell type mapping table for all clusters in the integrated dataset
all_integrated_cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/integration.allcluster2celltype.20200213.v3.tsv")
## input Alla's immune cell type assignment
immune_integrated_barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/Immune/30_aliquot_integration.barcode2celltype.20200330.AK.v1.tsv", data.table = F)
## input tumor subclustering cell type assignment
barcode2manualtumorcluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/annotate_barcode/annotate_barcode_with_manual_tumorsubcluster_id/20200324.v1/barcode2tumorsubclusterid.20200324.v1.tsv", data.table = F)

# map cell types for immune cell group ----------------------------------------------------------
## filter to immune cells only and select columns
immune_integrated_barcode2celltype_df <- immune_integrated_barcode2celltype_df %>%
  filter(!(Cell_type.detailed %in% c("Unknown", "Tumor", "Normal_Nephron_Epithelium", "Myofibroblasts", "Fibroblasts", "Endothelial cells"))) %>%
  mutate(manual_tumorsubcluster_id = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)
## see if the major cell types to correct
table(immune_integrated_barcode2celltype_df$Cell_type.shorter)
# B-cells              Basophils            CD4 T cells  CD4/CD8 proliferating  CD8 T cells activated  CD8 T cells exhausted 
# 223                     30                   3493                    468                    581                    996 
# cDC            Macrophages        Macrophages M2b Mixed myeloid/lymphoid               NK cells                    pDC 
# 454                  17669                     28                   1875                   1290                     47 
# Plasma                   Treg                    TRM 
# 217                    574                    518 

## correct B-cells using Alla's assignment
### see if there are more detailed cell type
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "B-cells"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "B-cells"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "B-cells"] <- "B-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "B-cells"] <- ""
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "B-cells"] <- ""

## correct Plasma cells using Alla's assignment
### see if there are more detailed cell type
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- "B-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- "Plasma cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- ""
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- "Plasma cells"

## correct Basophils using Alla's assignment
### see if there are more detailed cell type
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Basophils"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Basophils"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Basophils"] <- "Granulocytes"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Basophils"] <- "Basophils"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Basophils"] <- ""

## correct CD4+ T-cells using Alla's assignment
### see if there are more detailed cell type
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- "CD4+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- ""
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- "CD4+ T-cells"

## correct CD4/CD8 proliferating using Alla's assignment
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4/CD8 proliferating"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4/CD8 proliferating"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4/CD8 proliferating"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4/CD8 proliferating"] <- ""
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4/CD8 proliferating"] <- ""

## correct CD8 T cells activated using Alla's assignment
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- "CD8+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- ""
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- "CD8+ T-cells activated"

## correct CD8 T cells exhausted using Alla's assignment
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "CD8+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- ""
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "CD8+ T-cells exhausted"

## correct cDC exhausted using Alla's assignment
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "cDC"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "cDC"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "cDC"] <- "Monocytic lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "cDC"] <- "DC"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "cDC"] <- "cDC"

## correct pDC
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "pDC"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "pDC"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "pDC"] <- "Monocytic lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "pDC"] <- "DC"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "pDC"] <- "pDC"

## correct macrophage cell types
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages"] <- "Monocytic lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages"] <- "Macrophages"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages"] <- ""

## correct macrophage m2b cell types
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages M2b"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages M2b"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages M2b"] <- "Monocytic lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages M2b"] <- "Macrophages"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Macrophages M2b"] <- "Macrophages M2"

## correct mixed myeloid/lymphoid immune cells
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Mixed myeloid/lymphoid"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Mixed myeloid/lymphoid"] <- "Mixed myeloid/lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Mixed myeloid/lymphoid"] <- ""
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Mixed myeloid/lymphoid"] <- ""
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Mixed myeloid/lymphoid"] <- ""

## correct NK cells
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "NK cells"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "NK cells"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "NK cells"] <- "NK cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "NK cells"] <- ""
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "NK cells"] <- ""

## correct Tregs
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "CD4+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "Tregs"
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "Tregs"

## correct TRM
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "CD4+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "CD4+ memory T-cells"

# for tumor cells, label the unknown cells --------------------------------
## merge barcode to cluster info with cluster to cell type info
all_integrated_barcode2celltype_df <- merge(all_integrated_barcode2cluster_df, all_integrated_cluster2celltype_df,
                                            by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## format the barcode column and others to bind with other info later
all_integrated_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(Is_Normal_Nephron_Epithelium = ((ident == 17) & (Most_Enriched_Cell_Group == "Nephron_Epithelium"))) %>%
  mutate(integrated_barcode = barcode)
## get the integrated barcodes for tumor cells
tumor_barcode2celltype_df <- merge(all_integrated_barcode2celltype_df %>%
                                     select(orig.ident, individual_barcode, integrated_barcode,
                                            Most_Enriched_Cell_Group, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4), 
                                   barcode2manualtumorcluster_df %>%
                                     select(orig.ident, barcode, manual_cluster_id),
                                   by.x = c("orig.ident", "individual_barcode"),
                                   by.y = c("orig.ident", "barcode"),
                                   all.y = T)

## label the tumor cells
tumor_barcode2celltype_df$Most_Enriched_Cell_Group[!is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- "Nephron_Epithelium"
tumor_barcode2celltype_df$Most_Enriched_Cell_Type1[!is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- "Proximal tubule"
tumor_barcode2celltype_df$Most_Enriched_Cell_Type2[!is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""
tumor_barcode2celltype_df$Most_Enriched_Cell_Type3[!is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""
tumor_barcode2celltype_df$Most_Enriched_Cell_Type4[!is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""

## label the unknown cells within tumor cells
tumor_barcode2celltype_df$Most_Enriched_Cell_Group[is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- "Unknown"
tumor_barcode2celltype_df$Most_Enriched_Cell_Type1[is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""
tumor_barcode2celltype_df$Most_Enriched_Cell_Type2[is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""
tumor_barcode2celltype_df$Most_Enriched_Cell_Type3[is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""
tumor_barcode2celltype_df$Most_Enriched_Cell_Type4[is.na(tumor_barcode2celltype_df$manual_cluster_id)] <- ""

## add/remove columns
tumor_barcode2celltype_df <- tumor_barcode2celltype_df %>%
  mutate(is_malignant = T) %>%
  mutate(Cell_type.shorter = ifelse(Most_Enriched_Cell_Group == "Unknown", "Unknown", "Tumor cells")) %>%
  mutate(Cell_type.detailed = Cell_type.shorter) %>%
  mutate(Is_Normal_Nephron_Epithelium = F) %>%
  rename(manual_tumorsubcluster_id = manual_cluster_id) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# map stromal cells -------------------------------------------------------
stroma_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(Most_Enriched_Cell_Group == "Stroma") %>%
  mutate(Cell_type.shorter = Most_Enriched_Cell_Type1) %>%
  mutate(Cell_type.detailed = Most_Enriched_Cell_Type1)
## see what needs to be corrected
table(stroma_barcode2celltype_df$Most_Enriched_Cell_Type1)
table(stroma_barcode2celltype_df$Most_Enriched_Cell_Type2)
## make short cell type names for myofibroblasts
stroma_barcode2celltype_df$Cell_type.shorter[stroma_barcode2celltype_df$Most_Enriched_Cell_Type2 == "Myofibroblasts"] <- "Myofibroblasts"
stroma_barcode2celltype_df$Cell_type.detailed[stroma_barcode2celltype_df$Most_Enriched_Cell_Type2 == "Myofibroblasts"] <- "Myofibroblasts"
## add columns
stroma_barcode2celltype_df <- stroma_barcode2celltype_df %>%
  mutate(is_malignant = F) %>%
  mutate(manual_tumorsubcluster_id = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# map normal epithelial cells -------------------------------------------------------
normal_epithelial_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(Is_Normal_Nephron_Epithelium == T) %>%
  mutate(Cell_type.shorter = Most_Enriched_Cell_Type1) %>%
  mutate(Cell_type.detailed = Most_Enriched_Cell_Type2)
table(normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type1)
table(normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type2)
## add columns
normal_epithelial_barcode2celltype_df <- normal_epithelial_barcode2celltype_df %>%
  mutate(is_malignant = F) %>%
  mutate(manual_tumorsubcluster_id = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# map unknown cells -------------------------------------------------------
assigned_integrated_barcodes <- c(immune_integrated_barcode2celltype_df$integrated_barcode,
                                  tumor_barcode2celltype_df$integrated_barcode,
                                  stroma_barcode2celltype_df$integrated_barcode,
                                  normal_epithelial_barcode2celltype_df$integrated_barcode)
length(assigned_integrated_barcodes)
nrow(all_integrated_barcode2celltype_df)

unknown_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(!(integrated_barcode %in% assigned_integrated_barcodes)) %>%
  mutate(Cell_type.shorter = "Unknown") %>%
  mutate(Cell_type.detailed = "Unknown") %>%
  mutate(is_malignant = F) %>%
  mutate(manual_tumorsubcluster_id = NA)
unknown_barcode2celltype_df$Most_Enriched_Cell_Group <- "Unknown"
unknown_barcode2celltype_df$Most_Enriched_Cell_Type1 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type2 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type3 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type4 <- ""

unknown_barcode2celltype_df <- unknown_barcode2celltype_df %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# merge the immune and tumor cells and other info ------------------------------------
barcode2celltype_df <- rbind(immune_integrated_barcode2celltype_df,
                             tumor_barcode2celltype_df,
                             stroma_barcode2celltype_df,
                             normal_epithelial_barcode2celltype_df,
                             unknown_barcode2celltype_df)

# write output ------------------------------------------------------------
write.table(x = barcode2celltype_df, file = paste0(dir_out, "30_aliquot_integration.barcode2celltype.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
