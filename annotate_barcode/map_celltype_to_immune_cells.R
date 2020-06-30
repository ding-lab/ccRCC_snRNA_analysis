# Yige Wu @WashU Jun 2020

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
## input Alla's immune cell type assignment
immune_integrated_barcode2celltype_df <- fread(input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Immune/30_aliquot_integration.barcode2celltype.20200330.AK.v1.tsv", data.table = F)

# map cell types for immune cell group ----------------------------------------------------------
## filter to immune cells only and select columns
immune_integrated_barcode2celltype_df <- immune_integrated_barcode2celltype_df %>%
  filter(!(Cell_type.detailed %in% c("Unknown", "Tumor", "Normal_Nephron_Epithelium", "Myofibroblasts", "Fibroblasts", "Endothelial cells"))) %>%
  mutate(Id_TumorManualCluster = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)
nrow(immune_integrated_barcode2celltype_df)
## see if the major cell types to correct
table(immune_integrated_barcode2celltype_df$Cell_type.shorter)
# B-cells              Basophils            CD4 T cells  CD4/CD8 proliferating  CD8 T cells activated  CD8 T cells exhausted 
# 223                     30                   3493                    468                    581                    996 
# cDC            Macrophages        Macrophages M2b Mixed myeloid/lymphoid               NK cells                    pDC 
# 454                  17669                     28                   1875                   1290                     47 
# Plasma                   Treg                    TRM 
# 217                    574                    518 

## no unknown

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
immune_integrated_barcode2celltype_df$Cell_type.shorter[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Plasma"] <- "Plasma cells"


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
immune_integrated_barcode2celltype_df$Cell_type.shorter[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD4 T cells"] <- "CD4+ T-cells"

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
immune_integrated_barcode2celltype_df$Cell_type.shorter[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells activated"] <- "CD8+ T-cells activated"

## correct CD8 T cells exhausted using Alla's assignment
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "Lymphoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "CD8+ T-cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- ""
immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "CD8+ T-cells exhausted"
immune_integrated_barcode2celltype_df$Cell_type.shorter[immune_integrated_barcode2celltype_df$Cell_type.shorter == "CD8 T cells exhausted"] <- "CD8+ T-cells exhausted"

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

## correct TRM
table(immune_integrated_barcode2celltype_df$Cell_type.detailed[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"])
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "Myleoid lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "Monocytic lineage immune cells"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "Macrophages"
immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4[immune_integrated_barcode2celltype_df$Cell_type.shorter == "TRM"] <- "TRM"

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
immune_integrated_barcode2celltype_df$Cell_type.shorter[immune_integrated_barcode2celltype_df$Cell_type.shorter == "Treg"] <- "Tregs"



immune_integrated_barcode2celltype_df$Most_Enriched_Cell_Group <- "Immune"
# write output ------------------------------------------------------------
nrow(immune_integrated_barcode2celltype_df)
## [1] 33004
file2write <- paste0(dir_out, "Barcode2ImmuneCellType.", run_id, ".tsv")
write.table(x = immune_integrated_barcode2celltype_df, file = file2write, sep = '\t', quote = F, row.names = F)


