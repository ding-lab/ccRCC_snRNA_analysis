# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input clinical data
clinical_case_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data_snrna_samples/20211011.v1/snRNA_ccRCC_Clinicl_Table.20211011.v1.tsv")
## input FASTQ summary
fastq_summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.tsv")
fastq_summary_list <- readRDS(file = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.RSD")

# concatenate FASTQ summary -----------------------------------------------
fastq_summary_df <- NULL
for (file_tmp in names(fastq_summary_list)) {
  df_tmp <- fastq_summary_list[[file_tmp]]
  colname_lane <- colnames(df_tmp)[grepl(pattern = "Lane", x = colnames(df_tmp))]
  colname_samplename <- colnames(df_tmp)[grepl(pattern = "Sample", x = colnames(df_tmp))]
  
  df_tmp2 <- df_tmp[c("Flow Cell ID", "Index Sequence", "Completion Date", colname_samplename, colname_lane, "Library Name", "data_type", "file_dir")]
  colnames(df_tmp2) <- c("Flow Cell ID", "Index Sequence", "Completion Date", "Sample Name", "Lane Number", "Library Name", "data_type", "file_dir")
  df_tmp2[, "Completion Date"] <- as.character(df_tmp2[, "Completion Date"])
  fastq_summary_df <- rbind(fastq_summary_df, df_tmp2)
}

# preprocess fastq summary info -------------------------------------------
fastq_summary_ccRCC_df <-  fastq_summary_df %>%
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
  mutate(Location =  paste0("/diskmnt/primary/ccRCC_snRNA/", Aliquot.snRNA, "/")) %>%
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


final_df <- rbind(snRNA_metadata_df,
                  metadata_filtered_df %>%
                    filter(snATAC_used) %>%
                    mutate(library_strategy = "ATAC-seq") %>%
                    mutate(path_katmai = paste0("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Cell_Ranger/outputs/", Aliquot.snRNA, "/outs/possorted_bam.bam")) %>%
                    mutate(expected_cells = NA) %>%
                    rename(sample_type.sim = Sample_Type) %>%
                    mutate(snRNA_aliquot_id = NA) %>%
                    mutate(snATAC_aliquot_id = Aliquot.snRNA) %>%
                    mutate(snATAC_experiment_type = "snATAC") %>%
                    mutate(snRNA_experiment_type = NA) %>%
                    mutate(has_snATAC = "TRUE") %>%
                    mutate(has_snRNA = "FALSE") %>%
                    select(Case, snATAC_aliquot_id, snRNA_aliquot_id, snATAC_experiment_type, snRNA_experiment_type, has_snATAC, has_snRNA,
                           sample_name, sample_type.sim, library_strategy, expected_cells, path_katmai))
final_df <- final_df %>%
  rename(case_id = Case) %>%
  mutate(sample_type = ifelse(sample_type.sim == "Tumor", "Primary Tumor", "Normal adjacent tissue")) %>%
  mutate(Note = "NA") %>%
  mutate(Published = NA) %>%
  select(case_id, sample_type, snATAC_aliquot_id, snRNA_aliquot_id, snATAC_experiment_type, snRNA_experiment_type, has_snATAC, has_snRNA, Published)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GDC.ccRCC.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


