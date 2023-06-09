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
# metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
metadata_df <- fread(data.table = F, input = "~/Documents/Project/ccRCC_snRNA/Resouces/meta_data.20210809.v1.tsv")
## input FASTQ summary
fastq_summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.tsv")
fastq_summary_list <- readRDS("./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.RSD")
## input patient-level clinical data
patient_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data/20230310.v1/snRNA_ccRCC_Clinicl_Table.20230310.v1.tsv")
## input file sizes
filesizes_df <- fread(data.table = F, input = "~/Documents/Project/ccRCC_snRNA/Submission/CDS_upload/ccRCC.snATAC.fastq.gz.filesizes.txt", col.names = c("file_size", "path"))
## input md5
md5sum_df <- fread(data.table = F, input = "~/Documents/Project/ccRCC_snRNA/Submission/CDS_upload/extract_md5sum/20230316.v2/ccRCC.snRNA.snATAC.md5sum.20230316.v2.tsv")
md5sum_df <- fread(data.table = F, input = "~/Documents/Project/ccRCC_snRNA/Submission/CDS_upload/fastq.gz.md5sum.txt", col.names = c("md5sum", "file_path"), header = F)

# preprocess fastq summary info -------------------------------------------
# View(fastq_summary_list[[1]])
## this info is extracted from the sampleinfo.csv we grabbed recursively in the primary data folder
fastqinfo_summary_df <- NULL
for (idx in names(fastq_summary_list)) {
  fastqinfo_tmp_df <- fastq_summary_list[[idx]]
  if ("snATAC" %in% fastqinfo_tmp_df$data_type) { ## we are only extracting the snATAC FASTQ files
    colname_samplename <- colnames(fastqinfo_tmp_df)[grepl(pattern = "Sample", x = colnames(fastqinfo_tmp_df))]
    if ("File Name" %in% colnames(fastqinfo_tmp_df)) {
      fastqinfo_keep_tmp_df <- fastqinfo_tmp_df[, c("file_dir", "File Name", colname_samplename)]
    } else {
      colnames2process <- colnames(fastqinfo_tmp_df)[which(grepl(x = fastqinfo_tmp_df[1,], pattern = "fastq"))]
      fastqinfo_tmp_df1 <- fastqinfo_tmp_df[, c("file_dir", colname_samplename, colnames2process)]
      fastqinfo_tmp_df2 <- melt.data.table(data = data.table(fastqinfo_tmp_df1), measure.vars = colnames2process)
      fastqinfo_tmp_df2 <- as.data.frame(fastqinfo_tmp_df2)
      fastqinfo_keep_tmp_df <- fastqinfo_tmp_df2[, c("file_dir", "value", colname_samplename)]
    }
    colnames(fastqinfo_keep_tmp_df) <- c("file_dir", "File Name", "Sample Name")
    fastqinfo_summary_df <- rbind(fastqinfo_summary_df, fastqinfo_keep_tmp_df)
  }
}
unique(fastqinfo_summary_df$file_dir)
unique(fastqinfo_summary_df$`Sample Name`)
## so far the samples includes some samples not included in this study, we need to filter these out
fastq_summary_ccRCC_df <-  fastqinfo_summary_df %>%
  filter(file_dir != "2863682_new") %>%
  select(-file_dir) %>% # as the same sampleinfo.csv contains multiple samples, so I have copied the same file to different sample folders
  unique() %>%
  mutate(aliquot = str_split_fixed(string = `Sample Name`, pattern = "\\-", n = 3)[,2]) %>%
  mutate(aliquot = gsub(x = aliquot, pattern="N1", replacement=""))

# preprocess meta data --------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(snATAC_used) %>%
  select(Aliquot.snRNA.WU, Case, Sample, Sample_Type, Aliquot.snRNA)

# merge info-----------------------------------------------------
metadata_merged_df <- merge(x = metadata_filtered_df, 
                            y = fastq_summary_ccRCC_df, 
                            by.x = c("Aliquot.snRNA"), by.y = c("aliquot"), all.x = T)
metadata_add_df = metadata_merged_df[metadata_merged_df$Aliquot.snRNA == "CPT0000890002",]
metadata_add_df$`File Name` = c("TWFU-CPT0000890002-XBa1_S2_L002_I1_001.fastq.gz", "TWFU-CPT0000890002-XBa1_S2_L002_R1_001.fastq.gz")
metadata_merged_df = rbind(metadata_merged_df, metadata_add_df)
unique(metadata_merged_df$Aliquot.snRNA) ## 28 samples
table(metadata_merged_df$Aliquot.snRNA)
metadata_merged_df <- merge(x = metadata_merged_df, y = patient_clinical_df, by = c("Case"), all.x = T) ## this has to be added because gender info is required
filesizes_df$file_name <- sapply(filesizes_df$path, FUN = function(p) {
  p_split = str_split(p, pattern = "\\/")[[1]]
  filename = p_split[length(p_split)]
  return(filename)
})
filesizes_df = filesizes_df[!grepl(pattern = "2863682_new", x = filesizes_df$path),]
## add file sizes
metadata_merged_df <- merge(x = metadata_merged_df, 
                            y = filesizes_df,
                            by.x = c("File Name"), by.y = c("file_name"), all.x = T)
## add md5sum
md5sum_df$file_name.orig = sapply(X = md5sum_df$file_path, FUN = function(p) {
  p_split = str_split(string = p, pattern = "\\/")[[1]]
  filename_tmp = p_split[length(p_split)]
  return(filename_tmp)
})
md5sum_df = md5sum_df[!grepl(pattern = "2863682_new", x = md5sum_df$file_path),]

metadata_merged_df <- merge(x = metadata_merged_df, 
                            y = md5sum_df %>%
                              # mutate(file_name.orig = gsub(x = file_name, pattern = "\\.md5", replacement = "")) %>%
                              select(file_name.orig, md5sum, file_path),
                            by.x = c("File Name"), by.y = c("file_name.orig"), all.x = T)
## examine final numbers
unique(metadata_df$Case[metadata_df$snRNA_available & metadata_df$Case != "C3L-00359"])
unique(metadata_df$Sample[metadata_df$snRNA_available & metadata_df$Case != "C3L-00359"])

# create metadata fields --------------------------------------------------
cds_metadata_df <- metadata_merged_df %>%
  ## Study Top Level
  mutate(phs_accession="phs001287.v16.p6", 
         study_name="Washington University in St. Louis ccRCC snRNA-seq and snATAC-seq study",
         study_acronym="WashUccRCCsn",
         number_of_participants=25,
         number_of_samples=34,
         study_data_types="Genomic",
         experimental_strategy_and_data_subtype="RNA-seq",
         acl="['phs001287']",
         primary_investigator_name="Dr Li Ding",
         primary_investigator_email="lding@wustl.edu",
         co_primary_investigator_name="Dr Feng Chen",
         co_primary_investigator_email="fchen@wustl.edu",
         bioproject_accession="") %>%
  ## Funding information
  mutate(funding_agency="NCI",
         funding_source_program_name="",
         grant_id="R01HG009711;U24CA211006;U2CCA233303;U24CA210976",
         clinical_trial_system="",
         clinical_trial_identifier="",
         clinical_trial_arm="") %>%
  ## Additional Study Information (Study Summary)
  mutate(organism_species="Human",
         adult_or_childhood_study="adult",
         file_types_and_format="FastQ",
         size_of_data_being_uploaded="") %>%
  ## Participant (Clinical/Demographic - Patient Level)
  mutate(participant_id=Case,
         gender=Gender,
         race=ifelse(Race %in% c("White", "American indian or Alaska Native", "Black or African American", "Asian", "Native Hawaiian or Other Pacific islander"), Race, "Other"),
         ethnicity=ifelse(is.na(Ethnicity), 
                          "Unknown", 
                          ifelse(Ethnicity == "Not reported", 
                                 "Not Reported",
                                 ifelse(Ethnicity == "Not-Hispanic or Latino", 
                                        "Not Hispanic or Latino", 
                                        Ethnicity))),
         dbGaP_subject_id="",
         alternate_participant_id_1="",
         alternate_system_1="") %>%
  ## Sample
  mutate(sample_id=Case,
         sample_description="",
         biosample_accession="",
         sample_type=ifelse(Sample_Type == "Tumor", "Clear Cell Renal Cell Carcinoma - Kidney", "Normal Adjacent Tissue - Kidney")) %>%
  ## Additional Sample Information
  mutate(sample_tumor_status=ifelse(Sample_Type == "Tumor", "tumor", "normal"),
         sample_anatomic_site="Kidney",
         sample_age_at_collection="",
         derived_from_specimen=Sample) %>%
  ## File
  mutate(file_name=`File Name`,
         file_type="FASTQ",
         file_size=file_size, ## Required
         md5sum=md5sum,## Required
         file_url_in_cds="", ## Required after sending to CDS
         checksum_value="",
         checksum_algorithm="",
         file_mapping_level="Sample") %>%
  ## Sequencing-Specific Information
  mutate(library_id=`Sample Name`,
         library_strategy="ATAC-seq",
         library_source="GENOMIC SINGLE CELL",
         library_selection="",
         library_layout="paired-end",
         platform="ILLUMINA",
         instrument_model="",
         design_description="",
         reference_genome_assembly="GRCh38",
         custom_assembly_fasta_file_for_alignment="",
         bases="",
         number_of_reads="",
         coverage="",
         avg_read_length="",
         sequence_alignment_software="") %>%
  ## Diagnosis
  mutate(diagnosis_id="", # required but don't know what to enter 
         disease_type="Epithelial Neoplasm",
         primary_site="Kidney",
         vital_status="") %>%
  ## Additional Diagnosis Information
  mutate(age_at_diagnosis=Age,
         days_to_last_followup="",
         primary_diagnosis="Renal cell carcinoma, NOS",
         morphology="",
         tissue_or_organ_of_origin="Kidney",
         site_of_resection_or_biopsy="Kidney",
         tumor_grade="",
         tumor_stage_clinical_t=Primary_Tumor_Pathologic_Stage,
         tumor_stage_clinical_n="",
         tumor_stage_clinical_m="",
         progression_or_recurrence="",
         days_to_recurrence="",
         days_to_last_known_disease_status="",
         last_known_disease_status="") %>%
  select(phs_accession,
         study_name,
         study_acronym,
         number_of_participants,
         number_of_samples,
         study_data_types,
         experimental_strategy_and_data_subtype,
         acl,
         primary_investigator_name,
         primary_investigator_email,
         co_primary_investigator_name,
         co_primary_investigator_email,
         bioproject_accession,
         funding_agency,
         funding_source_program_name,
         grant_id,
         clinical_trial_system,
         clinical_trial_identifier,
         clinical_trial_arm,
         organism_species,
         adult_or_childhood_study,
         file_types_and_format,
         size_of_data_being_uploaded,
         participant_id,
         gender,
         ethnicity_race_ancestry_identified,
         race,
         ethnicity,
         dbGaP_subject_id,
         alternate_participant_id_1,
         alternate_system_1,
         sample_id,
         sample_description,
         biosample_accession,
         sample_type,
         sample_tumor_status,
         sample_anatomic_site,
         sample_age_at_collection,
         derived_from_specimen,
         file_name,
         file_type,
         file_size,
         md5sum,
         file_url_in_cds,
         checksum_value,
         checksum_algorithm,
         file_mapping_level,
         library_id,
         library_strategy,
         library_source,
         library_selection,
         library_layout,
         platform,
         instrument_model,
         design_description,
         reference_genome_assembly,
         custom_assembly_fasta_file_for_alignment,
         bases,
         number_of_reads,
         coverage,
         avg_read_length,
         sequence_alignment_software,
         diagnosis_id,
         disease_type,
         primary_site,
         vital_status,
         age_at_diagnosis,
         days_to_last_followup,
         primary_diagnosis,
         morphology,
         tissue_or_organ_of_origin,
         site_of_resection_or_biopsy,
         tumor_grade,
         tumor_stage_clinical_t,
         tumor_stage_clinical_n,
         tumor_stage_clinical_m,
         progression_or_recurrence,
         days_to_recurrence,
         days_to_last_known_disease_status,
         last_known_disease_status)

table(cds_metadata_df$race)
table(as.character(cds_metadata_df$ethnicity))
table(as.character(cds_metadata_df$sample_id))
cds_metadata_df = unique(cds_metadata_df)

sampleid_mapping_df = cds_metadata_df %>%
  rename(subject_id = participant_id) %>%
  select(subject_id, sample_id) %>%
  unique() %>%
  rename(SUBJECT_ID = subject_id) %>%
  rename(SAMPLE_ID = sample_id)

consent_df = sampleid_mapping_df %>%
  select(SUBJECT_ID) %>%
  unique() %>%
  mutate(CONSENT = 1)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CDS.snATAC.ccRCC.", run_id, ".tsv")
write.table(x = cds_metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "sample_mapping_for_dbGAP.", run_id, ".tsv")
write.table(x = sampleid_mapping_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "subject_consent_for_dbGAP.", run_id, ".tsv")
write.table(x = consent_df, file = file2write, quote = F, sep = "\t", row.names = F)

filecount_df = cds_metadata_df %>%
  select(sample_id) %>%
  table() %>%
  as.data.frame()
