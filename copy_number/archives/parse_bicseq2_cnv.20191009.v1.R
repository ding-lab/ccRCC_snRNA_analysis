# Yige Wu @WashU Oct 2019
## for parsing the ccRCC WGS-based CNV results from BIC-seq2


# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# cancatenate BIC-seq2 per patient cnv result ----------------------------------------------------
batchName <- "run_cases.CCRCC.MGI"
dir_tmp <- paste0("./Ding_Lab/Projects_Current/CPTAC/CPTACIII/cptac3_cnv/cptac3_cnv_shared_data/BICSEQ2/outputs/", batchName, "/annotation/")
cnv_file_names <- list.files(path = dir_tmp)
cnv_file_names <- cnv_file_names[grepl(pattern = "gene_level.log2.seg", x = cnv_file_names)]

bicseq2_cnv_tab <- NULL
for (cnv_file_name in cnv_file_names) {
  cnv_file_path <- paste0(dir_tmp, cnv_file_name)
  if (file.info(cnv_file_path)["size"] > 0) {
    bicseq2_cnv_tmp <- fread(input = cnv_file_path, data.table = F)
    case <- strsplit(x = cnv_file_name, split = "\\.")[[1]][1]
    
    bicseq2_cnv_tmp2bind <- data.frame(bicseq2_cnv_tmp[,"V5"])
    colnames(bicseq2_cnv_tmp2bind) <- case
    
    if (is.null(bicseq2_cnv_tab)) {
      bicseq2_cnv_tab <- data.frame(gene = bicseq2_cnv_tmp$V1)
    }
    
    bicseq2_cnv_tab <- cbind(bicseq2_cnv_tab, bicseq2_cnv_tmp2bind)
  }
}

write.table(x = bicseq2_cnv_tab, file = paste0(dir_out, "CPTAC3_ccRCC_Discovery_Set_BICSEQ2.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
