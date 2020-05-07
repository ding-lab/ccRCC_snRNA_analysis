# Yige Wu @WashU Feb 2020

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
## input HIF targets
hiftargets_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv", data.table = F)
## input EMT markers
markers_emt_df <- readxl::read_xlsx(path = "./Resources/Gene_Lists/EMT/EMT_Markers.xlsx", sheet = "Sheet1")
## input immune markers
markers_immune_df <- readxl::read_xlsx(path = "./Resources/Gene_Lists/Cytokines/Cytokine_Genes.20200504.v1.xlsx", sheet = "Sheet1")

# input genes by phenotype -------------------------------------------------------
markers_cancerstemcells_df <- data.frame(gene_symbol = c("CXCR4", "ENG", "PROM1", "CD44", "POU5F1", "NANOG", "ALDH1"),
                                         gene_function = NA,
                                         regulate_process = c("Chemostaxis", rep(NA, 6)),
                                         phenotype = "Stem cell identity")
unique(hiftargets_df$target_genesymbol)
markers_hypoxia_df <- data.frame(gene_symbol = c("HIF1A", "EPAS1", 
                                                 "HIF1A-AS2", 
                                                 hiftargets_df$target_genesymbol),
                                 gene_function = c(rep(x = "Transcription factor", 2), 
                                                      "lncRNA",
                                                      rep(NA, length(hiftargets_df$target_genesymbol))),
                                 regulate_process = c(rep("HIF Pathway", 2), NA, hiftargets_df$target_genefunction),
                                 phenotype = "Hypoxic response")
markers_hypoxia_df <- unique(markers_hypoxia_df)

# aggregate ---------------------------------------------------------------
colnames_merge <- colnames(markers_cancerstemcells_df)
markers_ith_df <- rbind(markers_cancerstemcells_df, 
                        markers_hypoxia_df, 
                        markers_emt_df[,colnames_merge],
                        markers_immune_df[,colnames_merge])

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "markergenes_by_intratumorheterogeneity_types.", run_id, ".tsv")
write.table(x = markers_ith_df, file = file2write, sep = "\t", quote = F, row.names = F)