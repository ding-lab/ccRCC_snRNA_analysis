# Yige Wu @WashU May 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(Signac)
library(rtracklayer)
# BiocManager::install("liftOver")
library(liftOver)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
dap_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/20210510/da_up_peaks_Tumor_vs_PT.annotated.20210510.tsv")
## input chain file
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

# filter to non-promoter and convert --------------------------------------
dap_nonpromoter_df <- dap_df %>%
  filter(Type != "Promoter") %>%
  filter(Count_up >= 5) %>%
  mutate(chr = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1])
grange_obj <- StringToGRanges(regions = dap_nonpromoter_df$peak)
grange_obj
seqlevelsStyle(grange_obj) = "UCSC" 
grange_hg19_obj = liftOver(x = grange_obj, chain = ch)
grange_hg19_obj %>% length()
num_ranges_vec <- sapply(grange_hg19_obj, function(x) {
  len <- length(x@ranges)
  return(len)
})
which(num_ranges_vec > 1)
which(num_ranges_vec == 0)
grange_hg19_obj[[1065]]
grange_hg19_obj = unlist(grange_hg19_obj)
grange_hg19_obj %>% length()
genome(grange_hg19_obj) = "hg19"
grange_hg19_obj
range_hg19_df <- grange_hg19_obj@ranges %>% as.data.frame()
colnames(range_hg19_df)[1:2] <- c("start.hg19", "end.hg19")

# combine with original file ----------------------------------------------
dap_nonpromoter_hg19_df <- dap_nonpromoter_df[rep(1:nrow(dap_nonpromoter_df), num_ranges_vec),]
dap_nonpromoter_hg19_df <- cbind(dap_nonpromoter_hg19_df, range_hg19_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Up_DAP.NonPromoter.withHg19.tsv")
write.table(x = dap_nonpromoter_hg19_df, file = file2write, quote = F, sep = "\t", row.names = F)
