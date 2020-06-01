# Yige Wu @WashU May 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(copynumber)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input the formatted GATK4 WXS CNV segments
segment_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/copy_number/format_gatk4_wxs_cnv_for_segment_plotting/20200511.v1/WXS_CNV_Somatic.Segments.20200511.v1.tsv")

# plot by each sample -----------------------------------------------------
## get aliquots to process
aliquots_all <- unique(segment_df$sampleID)
aliquots_all
## process each aliquot
for (aliquot_tmp in aliquots_all) {
  ## get the segment for this aliquot only
  plotdata_df <- segment_df %>%
    filter(sampleID == aliquot_tmp)
  
  chrs_plot <- unique(plotdata_df$chrom)
  for (chr_tmp in chrs_plot) {
    ## create output directory
    dir_out1 <- paste0(dir_out, "Chr", chr_tmp, "/")
    dir.create(dir_out1)
    
    ## plotChrom
    filepath_write <- paste0(dir_out1, aliquot_tmp, ".WXS_CNV_Somatic.plotChrom.", "Chr", chr_tmp, ".png")
    png(filename = filepath_write, width = 2000, height = 800, res = 150)
    plotChrom(segments = plotdata_df, ylab = "Log2(TumorCopy/NormalCopy)", seg.lwd = 0.7, ylim = c(-2, 2), chrom = chr_tmp, cyto.text = T)
    dev.off()
  }
}