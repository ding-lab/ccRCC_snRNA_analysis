# Yige Wu @WashU May 2021

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
## input motif mapped
motifs_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/Motifs_matched.Up_DAP.Overlap.Up_snRNA_bulk_DEGs.20210511.v1.tsv")
## input pathway annotation
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway__tumorcells_up_degs_msigdb_H_CP_sig_pathways/20210430.v1/DEG2Pathway.tsv")

# filter down to only motifs mapped to promoter/enhancer regions ----------
motifs_filtered_df <- motifs_df %>%
  rename(motif_region_type = Type) %>%
  rename(motif_mapped_gene = Gene) %>%
  filter(DAP_Type.strict == "Enhancer" | (DAP_Type.strict == "Promoter" & motif_region_type == "Promoter" & genesymbol_deg == motif_mapped_gene))

count_motifs_df <- motifs_filtered_df %>%
  select(motif, TF_name, motif_mapped_gene) %>%
  unique() %>%
  select(motif, TF_name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))

count_motifs_bydaptype_df <- motifs_filtered_df %>%
  select(motif, TF_name, motif_mapped_gene, DAP_Type.strict) %>%
  unique() %>%
  select(motif, TF_name, DAP_Type.strict) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))
count_motifs_promoter_df <- motifs_filtered_df %>%
  filter(DAP_Type.strict == "Promoter") %>%
  select(motif, TF_name, motif_mapped_gene) %>%
  unique() %>%
  select(motif, TF_name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))
count_motifs_enhancer_df <- motifs_filtered_df %>%
  filter(DAP_Type.strict == "Enhancer") %>%
  select(motif, TF_name, motif_mapped_gene) %>%
  unique() %>%
  select(motif, TF_name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))

# merge with pathway annotation and count ---------------------------------
motifs_annotated_df <- merge(x = motifs_filtered_df, y = gene2pathway_df, by.x = c("motif_mapped_gene"), by.y = c("GeneSymbol"), all.x = T)
count_motifs_bypathway_df <- motifs_annotated_df %>%
  select(motif, TF_name, GeneSet_Name, motif_mapped_gene) %>%
  unique() %>%
  select(motif, TF_name, GeneSet_Name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Motif_in_ProEnh_DAP_DEGs.", run_id, ".tsv")
write.table(x = motifs_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Count_Motif_in_DAP_DEGs.", run_id, ".tsv")
write.table(x = count_motifs_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Count_Motif_in_DAP_DEGs.ByDAPType.", run_id, ".tsv")
write.table(x = count_motifs_bydaptype_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Count_Motif_in_DAP_DEGs.ByPathway.", run_id, ".tsv")
write.table(x = count_motifs_bypathway_df, file = file2write, quote = F, sep = "\t", row.names = F)

