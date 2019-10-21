# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set run id ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set aliquot ID -----------------------------------------------------------
sample_id <- "CPT0086350004"

# load infercnv object------------------------------------------------------------
infercnv_obj = readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191005.v1/CPT0086350004.infercnv_obj.20191005.v1.RDS")
"VHL" %in% rownames(infercnv_obj@count.data)
"VHL" %in% rownames(infercnv_obj@expr.data)
"VHL" %in% rownames(infercnv_obj@gene_order)
summary(infercnv_obj@count.data["VHL", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["BAP1", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["SETD2", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["KDM5C", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["PTEN", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["MTOR", unlist(infercnv_obj@reference_grouped_cell_indices)])
summary(infercnv_obj@count.data["TP53", unlist(infercnv_obj@reference_grouped_cell_indices)])

cut_off_tmp <- 0.04
dir_out_cutoff <- paste0(dir_out, "cutoff", cut_off_tmp, "/")
dir.create(dir_out_cutoff)

sink(file = paste0(dir_out_cutoff, sample_id, ".inferCNV_Run_Cutoff", cut_off_tmp, ".", run_id, ".txt"))

infercnv_run_obj = infercnv::run(infercnv_obj,
                             cutoff=cut_off_tmp, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=dir_out_cutoff,
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=T)
sink()

infercnv_obj@reference_grouped_cell_indices
