# Yige Wu @WashU May 2020

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
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## set input directory
dir_cnv_files <- "./Resources/Bulk_Processed_Data/WXS_CNV_Somatic/b38/"

# input segment file by case ----------------------------------------------
## get the cases to process
cases_snrna <- unique(idmetadata_df$Case[idmetadata_df$Is_discovery_set == T & idmetadata_df$Sample_Type == "Tumor"])
cases_snrna
## get all the file names in the cnv directory
filenames_cnv <- list.files(dir_cnv_files)
## initiate super table
segment_df <- NULL
for (case_tmp in cases_snrna) {
  ## get the file name for this case
  filename_tmp <- filenames_cnv[grepl(pattern = case_tmp, x = filenames_cnv)]
  if (length(filename_tmp) == 0) {
    print(paste0(case_tmp, " has no CNV file!"))
  } else {
    ## get the path to the file
    filepath_tmp <- paste0(dir_cnv_files, filename_tmp)
    ## input cnv file
    cnv_tmp <- fread(data.table = F, input = filepath_tmp)
    ## get the aliquot.wu
    aliquot.wu_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Case == case_tmp & idmetadata_df$Is_discovery_set == T & idmetadata_df$Sample_Type == "Tumor"]
    ## format for plotting the cnv segments
    segment_tmp <- cnv_tmp %>%
      mutate(sampleID = aliquot.wu_tmp) %>%
      rename(chrom = Chromosome) %>%
      rename(start.pos = Start) %>%
      rename(end.pos = End) %>%
      rename(n.probes = Num_Probes) %>%
      mutate(mean = log2(Segment_Mean)) %>%
      mutate(arm = NA) %>%
      select(sampleID, chrom, arm, start.pos, end.pos, n.probes, mean)
    ## merge with the super table
    segment_df <- rbind(segment_df, segment_tmp)
  }
}
# write output ------------------------------------------------------------
filepath_write <- paste0(dir_out, "WXS_CNV_Somatic.Segments.", run_id, ".tsv")
write.table(x = segment_df, file = filepath_write, quote = F, sep = "\t", row.names = F)

