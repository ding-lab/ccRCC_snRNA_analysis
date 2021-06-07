# Yige Wu @WashU Jun 2021

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
cnv_gtf_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/WGS_CNV_Somatic/Michigan/c3-ccrcc-combined-cnvex-seg_v1.0.gtf", skip = 4)

# process -----------------------------------------------------------------
cnv_gtf_df$lr <- str_split_fixed(string = cnv_gtf_df$V9, pattern = "; lr", n = 2)[,2]
cnv_gtf_df$lr <- str_split_fixed(string = cnv_gtf_df$lr, pattern = ";", n = 2)[,1]
cnv_gtf_df$lr <- gsub(x = cnv_gtf_df$lr, pattern = '"', replacement = "")
cnv_gtf_df$lr <- as.numeric(as.vector(cnv_gtf_df$lr))

# write output ------------------------------------------------------------
filepath_write <- paste0(dir_out, "WXS_CNV_Somatic.Segments.", run_id, ".tsv")
write.table(x = segment_df, file = filepath_write, quote = F, sep = "\t", row.names = F)

