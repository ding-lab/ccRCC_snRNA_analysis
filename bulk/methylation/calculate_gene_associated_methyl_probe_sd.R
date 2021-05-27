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
## input methylation
methyl_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/DNA_methylation/CPTAC_ccRCC_discovery_tumor_methylation_betavalue_probe_level_v1.0.tsv")
## input methylation annotation
methyl_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")
## input clear-cell RCC classification
histology_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# preprocess --------------------------------------------------------------
## figure out which cases will be tested
caseids2test <- histology_df$CASE_ID[histology_df$Histologic_Type == "Clear cell renal cell carcinoma"]
length(caseids2test) ## 103

# filter methylation ------------------------------------------------------
nrow(methyl_df) ## 917437 rows
## some of probes are duplicated, and the values are the same across samples
methyl_df$Locus[!(methyl_df$Locus %in% methyl_anno_df$probeID)]
methyl_df$Locus[duplicated(methyl_df$Locus)]
# methyl_df[methyl_df$Locus == "cg00000714",]
## only use probes with gene annotated
probeids_w_gene <- methyl_anno_df$probeID[!is.na(methyl_anno_df$gene)]
length(probeids_w_gene) ## 706421
## filter
methyl_df <- methyl_df[methyl_df$Locus %in% probeids_w_gene,]
# methyl_df <- unique(methyl_df)
methyl_df <- methyl_df[!duplicated(methyl_df$Locus),]
## filter methylation matrix
colnames_methyl <- colnames(methyl_df)
colnames_methyl_new <- gsub(pattern = "\\-T", replacement = "", x = colnames_methyl)
colnames(methyl_df) <- colnames_methyl_new
methyl_df <- methyl_df[, c("Locus", caseids2test)]

# calculate ---------------------------------------------------------------
# test_df <- methyl_df[1:100000, -1]
# ## try the sapply
# start_time <- Sys.time()
# sd_test <- apply(test_df, 1, sd, na.rm = T)
# end_time <- Sys.time()
# end_time - start_time
# # Time difference of 3.219247 secs
sd_methyl <- apply(methyl_df[,-1], 1, sd, na.rm = T)
sd_methyl_df <- data.frame(probe_id = methyl_df$Locus, sd = sd_methyl, datapoints = rowSums(!is.na(methyl_df[,-1])))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Methylation_Probe.SD.", run_id, ".tsv")
write.table(x = sd_methyl_df, file = file2write, quote = F, sep = "\t", row.names = F)



  
  
  