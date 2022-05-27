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
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input ------------------------------------------------------
counts_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/Knockdown_Cell_Lines/all.gene_counts.tsv")

# create a DGEList object -------------------------------------------------
counts_mat <- as.matrix(counts_df[, colnames(counts_df)[grepl(pattern = "caki_1_cp|caki_1_control", x = colnames(counts_df))]])
dim(counts_mat)
# rownames(counts_mat) <- counts_df$ensembl_gene_id
dgList <- DGEList(counts=counts_mat, genes=counts_df[,1:7])

# Normalization -----------------------------------------------------------
countsPerMillion <- cpm(dgList); head(countsPerMillion); summary(countsPerMillion)
countCheck <- countsPerMillion > 1
dgList <- dgList[rowSums(countCheck) >= 1,]
dgList <- calcNormFactors(dgList, method="TMM")
countsPerMillion <- cpm(dgList);

# Setting up the Model ----------------------------------------------------
# sampleGroup <- mapvalues(x = colnames(dgList), 
#                          from = c("sample.caki_1_mxi1_c1_e1","sample.caki_1_pcsk6_c1_e1",
#                                                         "sample.caki_1_cp_c2_e1","sample.caki_1_cp_c1_e1",
#                                                         "sample.caki_1_control_e1",
#                                                         "sample.rcc4_klf9_c3_e1","sample.rcc4_klf9_c2_e1",
#                                                         "sample.rcc4_mxi1_c2_e1","sample.rcc4_mxi1_c1_e1",
#                                                         "sample.rcc4_scrambled_e1", "sample.rcc4_control_e1"), 
#                          to = c("caki1_mxi1", "caki1_pcsk6", "caki1_cp", "caki1_cp", "caki1_control", "rcc4_klf9", "rcc4_klf9", "rcc4_mxi1", "rcc4_mxi1", "rcc4_scrambled", "rcc4_control"))
# sampleGroup <- factor(x = sampleGroup, levels = c("caki1_cp","caki1_control", "caki1_mxi1","caki1_pcsk6","rcc4_klf9","rcc4_mxi1","rcc4_scrambled","rcc4_control"))
sampleGroup <- mapvalues(x = colnames(dgList), 
                         from = c("sample.caki_1_cp_c2_e1","sample.caki_1_cp_c1_e1",
                                  "sample.caki_1_control_e1"), 
                         to = c("caki1_cp", "caki1_cp", "caki1_control"))
sampleGroup <- factor(x = sampleGroup, levels = c("caki1_cp","caki1_control"))

designMat <- model.matrix(~sampleGroup)
designMat

#  Estimating Dispersions -------------------------------------------------
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
file2write <- paste0(dir_out, "BCV.pdf")
pdf(file2write, width = 15, height = 15)
plotBCV(dgList)
dev.off()

# 2.8 Differential Expression ---------------------------------------------
fit <- glmQLFit(dgList, designMat)

qlf <- glmQLFTest(fit, contrast = c(1, -1))
result_df <- topTags(object = qlf, n = nrow(qlf$table))$table
nrow(result_df)
## count # of DEGs
result_df %>%
  filter(logFC < -1 & FDR < 0.05) %>%
  nrow()
result_df %>%
  filter(logFC > 1 & FDR < 0.05) %>%
  nrow()
result_df %>%
  filter(logFC < -10 & FDR < 0.05) %>%
  nrow()

result_df %>%
  filter(FDR > 0.05) %>%
  nrow()

## write output
file2write <- paste0(dir_out, "Caki1_CP_vs_Caki1_NT.DEGs.", run_id, ".tsv")
write.table(x = result_df, file = file2write, row.names = F, quote = F, sep = "\t")

