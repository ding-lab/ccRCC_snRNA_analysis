# Yige Wu @WashU May 2022
## reference: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf
## how they deal with count of 0 to log2CPM: https://support.bioconductor.org/p/107719/

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "edgeR"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}

# input ------------------------------------------------------
counts_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/other/preprocess/intersect_shRNA_lines_w_NT_lines_RNAseq/20220603.v1/Cell_Lines.gene_counts.20220603.v1.tsv")

# create a DGEList object -------------------------------------------------
# counts_mat <- as.matrix(counts_df[, colnames(counts_df)[grepl(pattern = "sample", x = colnames(counts_df))]])
counts_mat <- as.matrix(counts_df[, c("sample.caki_1_cp_c2_e1", "sample.caki_1_cp_c1_e1", "sample.caki_1_control_e1", "sample.rcc4_klf9_c2_e1", "sample.rcc4_mxi1_c2_e1",
                                      "sample.rcc4_scrambled_e1", "sample.dr_caki_1_rna", "sample.skrc42+bap1_e1", "sample.skrc42+emptyvector_e1")])
# rownames(counts_mat) <- counts_df$ensembl_gene_id
dgList <- DGEList(counts=counts_mat, genes=counts_df[,colnames(counts_df)[!grepl(pattern = "sample", x = colnames(counts_df))]])

# Normalization -----------------------------------------------------------
countsPerMillion <- cpm(dgList); head(countsPerMillion)
dgList <- calcNormFactors(dgList, method="TMM")
countsPerMillion_tmm <- cpm(dgList); head(countsPerMillion_tmm)
countsPerMillion_tmm_log2 <- cpm(dgList, log = T); head(countsPerMillion_tmm_log2)


# add information ---------------------------------------------------------
countsPerMillion_tmm_df <- cbind(counts_df[,1:7], data.frame(countsPerMillion_tmm))
countsPerMillion_tmm_log2_df <- cbind(counts_df[,1:7], data.frame(countsPerMillion_tmm_log2))


# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## write output
file2write <- paste0(dir_out, "CPM.TMM_normalized.Manuscript_Cell_Lines.", run_id, ".tsv")
write.table(x = countsPerMillion_tmm_df, file = file2write, row.names = F, quote = F, sep = "\t")
file2write <- paste0(dir_out, "Log2CPM.TMM_normalized.All_Cell_Lines.", run_id, ".tsv")
write.table(x = countsPerMillion_tmm_log2_df, file = file2write, row.names = T, quote = F, sep = "\t")

