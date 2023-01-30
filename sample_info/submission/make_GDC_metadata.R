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
clinical_case_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data_snrna_samples/20211011.v1/snRNA_ccRCC_Clinicl_Table.20211011.v1.tsv")
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
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, Sample, snRNA_available, snATAC_used, Sample_Type, Aliquot.snRNA)

# process snRNA data-----------------------------------------------------------------
snRNA_metadata_df <- merge(x = metadata_filtered_df, 
                            y = fastq_summary_ccRCC_df %>%
                              filter(data_type == "snRNA"), 
                            by.x = c("Aliquot.snRNA"), by.y = c("aliquot"), all.x = T)

snRNA_metadata_df <- snRNA_metadata_df %>%
  rename(case_id = Case) %>%
  rename(specimen_id = Sample) %>%
  rename(sample_type.sim = Sample_Type) %>%
  mutate(sample_type = ifelse(sample_type.sim == "Tumor", "Primary Tumor", "Normal adjacent tissue")) %>%
  mutate(Note = sample_type.sim) %>%
  mutate(snRNA_aliquot_id = Aliquot.snRNA) %>%
  mutate(snATAC_aliquot_id = NA) %>%
  mutate(snATAC_experiment_type = NA) %>%
  mutate(snRNA_experiment_type = "snRNA") %>%
  mutate(has_snATAC = "FALSE") %>%
  mutate(has_snRNA = "TRUE") %>%
  mutate(Published = "") %>%
  mutate(Status = "Ready for upload to GDC") %>%
  mutate(Location =  paste0("/diskmnt/primary/ccRCC_snRNA/", Aliquot.snRNA, "/BAM/")) %>%
  mutate(Path_BAM = paste0("/diskmnt/primary/ccRCC_snRNA/", Aliquot.snRNA, "/BAM/possorted_genome_bam.bam")) %>%
  rename(`Sequencing lane` = `Lane Number`) %>%
  rename(`Sample name` = `Sample Name`) %>%
  mutate(expected_cells = 6000) %>%
  mutate(`Sequence read length` = 150) %>%
  rename(`Number of expected cells` = expected_cells) %>%
  mutate(Sequencing_data_text = str_split_fixed(string = `Completion Date`, pattern = " ", n = 2)[,1]) %>%
  mutate(year = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                       str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,3],
                       substr(x = str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,1], start = 3, stop = 4))) %>%
  mutate(month = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                       str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,1],
                       str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,2])) %>%
  mutate(day = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                        str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,2],
                        str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,3])) %>%
  mutate(`Sequencing date` = paste0(as.numeric(month), "/", day, "/", year)) %>%
  mutate(Protocol = "10x_SC_3'_GEX_V3.1") %>%
  select(case_id, specimen_id, sample_type, Note,
         snATAC_aliquot_id, snRNA_aliquot_id, snATAC_experiment_type, snRNA_experiment_type, has_snATAC, has_snRNA,
         Published, Status, Location, Path_BAM, `Flow Cell ID`, `Index Sequence`, 
         `Sequencing lane`, `Sequence read length`, `Sample name`,
         `Number of expected cells`, `Sequencing date`, Protocol)

# process snATAC data -----------------------------------------------------
snATAC_metadata_df <- merge(x = metadata_filtered_df %>%
                              filter(snATAC_used), 
                            y = fastq_summary_ccRCC_df %>%
                              filter(data_type == "snATAC"), 
                            by.x = c("Aliquot.snRNA"), by.y = c("aliquot"), all.x = T)

snATAC_metadata_df <- snATAC_metadata_df %>%
  rename(case_id = Case) %>%
  rename(specimen_id = Sample) %>%
  rename(sample_type.sim = Sample_Type) %>%
  mutate(sample_type = ifelse(sample_type.sim == "Tumor", "Primary Tumor", "Normal adjacent tissue")) %>%
  mutate(Note = sample_type.sim) %>%
  mutate(snATAC_aliquot_id = Aliquot.snRNA) %>%
  mutate(snRNA_aliquot_id = NA) %>%
  mutate(snATAC_experiment_type = "snATAC") %>%
  mutate(snRNA_experiment_type = NA) %>%
  mutate(has_snRNA = "FALSE") %>%
  mutate(has_snATAC = "TRUE") %>%
  mutate(Published = "") %>%
  mutate(Status = "Ready for upload to GDC") %>%
  mutate(Location =  paste0("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Cell_Ranger/outputs/", Aliquot.snRNA, "/outs/")) %>%
  mutate(Path_BAM = paste0("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Cell_Ranger/outputs/", Aliquot.snRNA, "/outs/possorted_bam.bam")) %>%
  rename(`Sequencing lane` = `Lane Number`) %>%
  rename(`Sample name` = `Sample Name`) %>%
  mutate(expected_cells = 6000) %>%
  mutate(`Sequence read length` = 150) %>%
  rename(`Number of expected cells` = expected_cells) %>%
  mutate(Sequencing_data_text = str_split_fixed(string = `Completion Date`, pattern = " ", n = 2)[,1]) %>%
  mutate(year = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                       str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,3],
                       substr(x = str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,1], start = 3, stop = 4))) %>%
  mutate(month = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                        str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,1],
                        str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,2])) %>%
  mutate(day = ifelse(grepl(pattern = "\\/", x = Sequencing_data_text), 
                      str_split_fixed(string = Sequencing_data_text, pattern = "\\/", n = 3)[,2],
                      str_split_fixed(string = Sequencing_data_text, pattern = "\\-", n = 3)[,3])) %>%
  mutate(`Sequencing date` = paste0(as.numeric(month), "/", day, "/", year)) %>%
  mutate(Protocol = "10x_SC_ATAC_SEQ") %>%
  select(case_id, specimen_id, sample_type, Note,
         snATAC_aliquot_id, snRNA_aliquot_id, snATAC_experiment_type, snRNA_experiment_type, has_snATAC, has_snRNA,
         Published, Status, Location, Path_BAM, `Flow Cell ID`, `Index Sequence`, 
         `Sequencing lane`, `Sequence read length`, `Sample name`,
         `Number of expected cells`, `Sequencing date`, Protocol)

# combine -----------------------------------------------------------------
final_df <- rbind(snRNA_metadata_df, snATAC_metadata_df)
final_df$`Sequencing date`[final_df$`Sequencing date` == "NA//"] <- NA

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GDC.ccRCC.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


