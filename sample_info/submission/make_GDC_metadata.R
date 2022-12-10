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
fastq_summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221210.v1/ccRCC.snRNA.snATAC.fastq.summary.20221210.v1.tsv")

# preprocess fastq summary info -------------------------------------------
fastq_summary_df <-  fastq_summary_df %>%
  unique() %>%
  rename(Library_Name = `Library Name`) %>%
  mutate(aliquot = str_split_fixed(string = Library_Name, pattern = "\\-", n = 3)[,2])


# preprocess meta data --------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, snRNA_available, snATAC_used, Sample_Type, Aliquot.snRNA)

# process snRNA data-----------------------------------------------------------------
snRNA_metadata_df <- metadata_filtered_df %>%
  mutate(path_katmai = paste0("/diskmnt/primary/ccRCC_snRNA/", Aliquot.snRNA, "/BAM/possorted_genome_bam.bam")) %>%
  mutate(expected_cells = 6000) %>%
  rename(sample_type.sim = Sample_Type) %>%
  mutate(snRNA_aliquot_id = Aliquot.snRNA) %>%
  mutate(snATAC_aliquot_id = NA) %>%
  mutate(snATAC_experiment_type = NA) %>%
  mutate(snRNA_experiment_type = "snRNA") %>%
  mutate(has_snATAC = "FALSE") %>%
  mutate(has_snRNA = "TRUE") %>%
  select(Case, snATAC_aliquot_id, snRNA_aliquot_id, snATAC_experiment_type, snRNA_experiment_type, has_snATAC, has_snRNA,
         sample_type.sim, expected_cells, path_katmai)

snRNA_metadata_df2 <- merge(x = snRNA_metadata_df, y = fastq_summary_df %>%
                              filter(data_type == "snRNA"), by.x = c("snRNA_aliquot_id"), by.y = c("aliquot"), all.x = T)

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


