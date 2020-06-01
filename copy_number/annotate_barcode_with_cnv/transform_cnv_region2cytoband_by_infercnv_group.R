# Yige Wu @WashU May 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library("D3GB")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the CNV regions
group2cnv_region_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/unite_infercnv_group2cnv_region/20200511.v1/InferCNV_Group2CNV_Region.20200511.v1.tsv")
## edit Cytoband
GRCh38.bands <- GRCh38.bands %>%
  mutate(cytoband_name = paste0(chr, name))
## make function
chr_region2cytoband <- function(chr_vec, position_vec) {
  if (length(chr_vec) == length(position_vec)) {
    
    bandname_vec <- NULL
    for (i in 1:length(chr_vec)) {
      chr_tmp <- chr_vec[i]
      chr_tmp <- gsub(pattern = "chr", x = chr_tmp, replacement = "")
      pos_tmp <- position_vec[i]
      bandname_tmp <- GRCh38.bands[(GRCh38.bands$chr == chr_tmp) & (GRCh38.bands$start <= pos_tmp) & (GRCh38.bands$end > pos_tmp), "cytoband_name"]
      ## in case there is none for the cytoband
      bandname_tmp <- ifelse(length(bandname_tmp) == 1, as.vector(bandname_tmp), NA)
      bandname_vec <- c(bandname_vec, bandname_tmp)
    }
    return(bandname_vec)
  } else {
    return("Error:the length of chromosome number and position is not the same!")
  }
}


## make function fo mapping chromsomal position to cytoband
nrow(group2cnv_region_df)
group2cnv_region_df$start_cytoband <- chr_region2cytoband(group2cnv_region_df$chr, group2cnv_region_df$start)
group2cnv_region_df$end_cytoband <- chr_region2cytoband(group2cnv_region_df$chr, group2cnv_region_df$end)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CNV_Region2Cytoband_by_InferCNV_Group", run_id, ".tsv")
write.table(x = group2cnv_region_df, file = file2write, sep = "\t", row.names = F, quote = F)

