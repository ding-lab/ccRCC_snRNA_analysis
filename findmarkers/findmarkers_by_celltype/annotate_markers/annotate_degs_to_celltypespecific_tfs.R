# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type specific TF motifs
motif2cellgroup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_celltype_enriched_motifs_for_snRNA_degs/20200908.v2/CellType2TF.Top100byCellGroup.Enriched_Motifs.chromvar.MergedObj.byCell_type.20200908.v2.tsv")
## specify the top n degs
n_top <- 200
## input the DEGs
genes_df <- fread(input = paste0("./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200908.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.Top", n_top, "avg_logFC.tsv"), data.table = F)
## input TF-target table
target2tf_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/Transcriptional/omnipathdb.transcriptional.20200908.txt")


# get all the TFs for the cell type specific DEGs -------------------------
gene2tf_df <- merge(x = genes_df, y = target2tf_df, by.x = c("gene"), by.y = c("target_genesymbol"), all.x = T)
gene2tf_df <- gene2tf_df %>%
  arrange(cluster, gene)
gene2tf_knowntf_df <- gene2tf_df %>%
  filter(!is.na(source_genesymbol)) 


# merge with cell type specific TFs ---------------------------------------
## merge with gene2tf table
gene2motifenrichedtf_df <- merge(x = gene2tf_knowntf_df, y = motif2cellgroup_df, by.x = c("source_genesymbol", "cluster"), by.y = c("tf_genesymbol", "Cell_group.detailed"))
## filter by interaction direction
gene2motifenrichedtf_df <- gene2motifenrichedtf_df %>%
  filter(is_inhibition == 0)

# annotate gene to different levels of evidence ---------------------------
gene_anno_df <- data.frame(target_genesymbol = unique(genes_df$gene))
gene_anno_df$TF_evidence_level <- ifelse(gene_anno_df$target_genesymbol %in% gene2motifenrichedtf_df$gene,
                                         ifelse(gene_anno_df$target_genesymbol %in% gene2motifenrichedtf_df$gene[gene2motifenrichedtf_df$is_stimulation == 1],
                                                "Stimulated Target of Cell-Type-Motif-Enriched TFs",
                                                "Target of Cell-Type-Motif-Enriched TFs; Direction unknown"),
                                         "Not a Target")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellTypeDEGTop", n_top, "2CellTypeMotifEnrichedTFs.", run_id, ".tsv")
write.table(x = gene2motifenrichedtf_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "CellTypeDEGTop", n_top, "2CellTypeMotifEnrichedTFs.", "Annotation.", run_id, ".tsv")
write.table(x = gene_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)

