# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input clinical data
clinical_case_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data_snrna_samples/20211011.v1/snRNA_ccRCC_Clinicl_Table.20211011.v1.tsvsnRNA_ccRCC_Clinicl_Table.20211011.v1.tsv")
## input FASTQ summary
fastq_summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.tsv")

# preprocess fastq summary info -------------------------------------------
fastq_summary_ccRCC_df <-  fastq_summary_df %>%
  filter(file_dir != "2863682_new") %>%
  select(-file_dir) %>%
  unique() %>%
  rename(Library_Name = `Library Name`) %>%
  mutate(aliquot = str_split_fixed(string = `Sample Name`, pattern = "\\-", n = 3)[,2]) %>%
  mutate(aliquot = ifelse(aliquot == "C3L", str_split_fixed(string = `Sample Name`, pattern = "\\-", n = 5)[,4], aliquot))
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$Library_Name == "TWCE-C3L-00079-01-CPT000126_0013-lib1"] <- "CPT0001260013"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$Library_Name == "TWCE-C3L-00908-01-CPT008635_0004-lib1"] <- "CPT0086350004"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT007513004"] <- "CPT0075130004"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$Library_Name == "TWCE-CPT007514-0002-lib1"] <- "CPT0075140002"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$Library_Name == "TWCE-C3N-00733-CPT0025890002-lib1"] <- "CPT0025890002"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$Library_Name == "TWCE-C3N-00733-CPT0010110013-lib1"] <- "CPT0010110013"
# CPT0001220012 snRNA FASTQ is downloaded at /storage1/fs1/m.wyczalkowski/Active/Primary/HTAN.backup/xfer.genome.wustl.edu/gxfer1/50599374701002, but the Samplemap.csv file does not have info on this sample
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0000870003N1"] <- "CPT0000870003"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0010100001N1"] <- "CPT0010100001"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0023690004N1"] <- "CPT0023690004"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0025110004N1"] <- "CPT0025110004"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0075170013N1"] <- "CPT0075170013"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0079410004N1"] <- "CPT0079410004"
fastq_summary_ccRCC_df$aliquot[fastq_summary_ccRCC_df$aliquot == "CPT0086820004N1"] <- "CPT0086820004"
fastq_summary_ccRCC_df$`Flow Cell ID`[fastq_summary_ccRCC_df$aliquot %in% c("CPT0001260013", "CPT0000880001")] <- "HNWJHDRXX,HNW3FDRXX"
fastq_summary_ccRCC_df$`Completion Date`[fastq_summary_ccRCC_df$aliquot %in% c("CPT0001260013", "CPT0000880001")] <- "10/6/20"
fastq_summary_ccRCC_df$`Completion Date`[!grepl(pattern = "\\-|\\/", x = fastq_summary_ccRCC_df$`Completion Date`)] <- NA
fastq_summary_ccRCC_df <- unique(fastq_summary_ccRCC_df)

# preprocess meta data --------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(snRNA_available) %>%
  filter(Case != "C3L-00359") %>%
  select(Aliquot.snRNA.WU, Case, Sample, snRNA_available, snATAC_used, Sample_Type, Aliquot.snRNA)

# process snRNA data-----------------------------------------------------------------
snRNA_metadata_df <- metadata_filtered_df %>%
  mutate(library_name = paste0(Aliquot.snRNA, "_snRNA"),
         title = "Kidney, snRNAseq, CPTAC",
         organism = "Human",
         cell_type = "tissue",
         tissue = "Kidney",
         molecular = "nuclear RNA",
         single_or_paired_end = "paired-end",
         instrument_model = "Illumina HiSeq 4000/NovaSeq 6000",
         description = "10x Genomics",
         processed_data_file.1 = paste0(Aliquot.snRNA, "/outs/raw_feature_bc_matrix/features.tsv.gz"),
         processed_data_file.2 = paste0(Aliquot.snRNA, "/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),
         processed_data_file.3 = paste0(Aliquot.snRNA, "/outs/raw_feature_bc_matrix/matrix.mtx.gz"),
         raw_file = "Will be submitted to dbGaP") %>%
  select(library_name,
         title,
         organism,
         cell_type,
         tissue,
         molecular,
         single_or_paired_end,
         instrument_model,
         description,
         processed_data_file.1,
         processed_data_file.2,
         processed_data_file.3,
         raw_file)


# process snATAC data -----------------------------------------------------
snATAC_metadata_df <- metadata_filtered_df %>%
  filter(snATAC_used) %>%
  mutate(library_name = paste0(Aliquot.snRNA, "_snATAC"),
         title = "Kidney, snATACseq, CPTAC",
         organism = "Human",
         cell_type = "tissue",
         tissue = "Kidney",
         molecular = "genomic DNA",
         single_or_paired_end = "paired-end",
         instrument_model = "Illumina HiSeq 4000/NovaSeq 6000",
         description = "10x Genomics",
         processed_data_file.1 = paste0(Aliquot.snRNA, "/outs/raw_peak_bc_matrix/peaks.bed"),
         processed_data_file.2 = paste0(Aliquot.snRNA, "/outs/raw_peak_bc_matrix/barcodes.tsv"),
         processed_data_file.3 = paste0(Aliquot.snRNA, "/outs/raw_peak_bc_matrix/matrix.mtx"),
         raw_file = "Will be submitted to CDS") %>%
  select(library_name,
         title,
         organism,
         cell_type,
         tissue,
         molecular,
         single_or_paired_end,
         instrument_model,
         description,
         processed_data_file.1,
         processed_data_file.2,
         processed_data_file.3,
         raw_file)

# combine -----------------------------------------------------------------
final_df <- rbind(snRNA_metadata_df, snATAC_metadata_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GEO.ccRCC.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


