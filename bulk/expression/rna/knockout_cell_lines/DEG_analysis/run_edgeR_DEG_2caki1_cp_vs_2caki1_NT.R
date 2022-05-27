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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input ------------------------------------------------------
counts_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/preprocess/combine_shRNA_lines_w_NT_lines/20220517.v1/Cell_Lines.gene_counts.20220517.v1.tsv")

# create a DGEList object -------------------------------------------------
counts_mat <- as.matrix(counts_df[, colnames(counts_df)[grepl(pattern = "caki_1_cp|caki_1_control|caki_1_rna", x = colnames(counts_df))]])
dim(counts_mat)
# rownames(counts_mat) <- counts_df$ensembl_gene_id
dgList <- DGEList(counts=counts_mat, genes=counts_df[,1:9])

# Normalization -----------------------------------------------------------
countsPerMillion <- cpm(dgList); head(countsPerMillion); summary(countsPerMillion)
countCheck <- countsPerMillion > 1
"CP" %in% counts_df$external_gene_name[rowSums(countCheck) >= 1]
"CP" %in% counts_df$external_gene_name[rowSums(countCheck) >= 2]
"VEGFA" %in% counts_df$external_gene_name[rowSums(countCheck) >= 2]
dgList <- dgList[rowSums(countCheck) >= 2,]
dgList <- calcNormFactors(dgList, method="TMM")

# Setting up the Model ----------------------------------------------------
colnames(dgList)
sampleGroup <- mapvalues(x = colnames(dgList), 
                         from = c("sample.caki_1_cp_c2_e1","sample.caki_1_cp_c1_e1",
                                  "sample.caki_1_control_e1", "sample.dr_caki_1_rna"), 
                         to = c("caki1_cp", "caki1_cp", "caki1_nt", "caki1_nt"))
sampleGroup <- factor(x = sampleGroup, levels = c("caki1_nt","caki1_cp"))
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
qlf <- glmQLFTest(fit) ## caki1_cp vs. caki1_nt
result_df <- topTags(object = qlf, n = nrow(qlf$table))$table
nrow(result_df) #3 13945
## count # of DEGs
result_df %>%
  filter(logFC < 0 & FDR < 0.05) %>%
  nrow() ## 2286
result_df %>%
  filter(logFC > 0 & FDR < 0.05) %>%
  nrow() ## 2114
result_df %>%
  filter(FDR > 0.05) %>%
  nrow()
result_df %>%
  filter(FDR < 0.05) %>%
  nrow()
result_df %>%
  filter(logFC < 0 & FDR < 0.05) %>%
  arrange(logFC) %>%
  View()

## write output
file2write <- paste0(dir_out, "Caki1_CP_vs_Caki1_NT.DEGs.", run_id, ".tsv")
write.table(x = result_df, file = file2write, row.names = F, quote = F, sep = "\t")
