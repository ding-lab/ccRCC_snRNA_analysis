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
## input FASTQ info
md5_df = fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_filepath_md5sum/20230609.v2/ccRCC.snRNA.snATAC.md5sum.20230609.v2.tsv")
fastq_paths_df = fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_filepath_md5sum/20230609.v2/ccRCC.snRNA.snATAC.FASTQ.paths.20230609.v2.tsv")
fastq_info_samplemap_df = fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_filepath_md5sum/20230609.v2/ccRCC.snRNA.snATAC.FASTQ.detailedinfo.20230609.v2.tsv")

# preprocess meta data --------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, Sample, snRNA_available, snATAC_used, Sample_Type, Aliquot.snRNA)
aliquots_all = unique(metadata_filtered_df$Aliquot.snRNA)

# preprocess fastq summary info -------------------------------------------
## process MD5
md5_df = md5_df %>%
  rename(file_name.md5 = file_name) %>%
  mutate(file_name = gsub(x = file_name.md5, pattern = "\\.md5", replacement = "")) %>%
  rename(file_path.md5 = file_path) %>%
  mutate(file_path = gsub(x = file_path.md5, pattern = "\\.md5", replacement = "")) %>%
  mutate(file_path = gsub(x = file_path, pattern = "\\/MD5\\/", replacement = "/FASTQ/"))
## merge data 1
fastq_paths_df$file_name = sapply(fastq_paths_df$file_path, function(x) {
  str_vec = str_split(x, pattern = "\\/")[[1]]
  file_name = str_vec[length(str_vec)]
  return(file_name)
})
fastq_paths_df = data.frame(fastq_paths_df)
fastq_info_df = merge(x = fastq_paths_df, y = md5_df, by = c("file_name", "file_path"), all.x = T)

fastq_info_df = fastq_info_df %>%
  filter(grepl(x = file_path, pattern = "ccRCC_snRNA")) %>%
  mutate(aliquot = str_split_fixed(file_path, "\\/", 6)[,5]) %>%
  filter(aliquot %in% metadata_filtered_df$Aliquot.snRNA)
length(unique(fastq_info_df$aliquot))
aliquots_all[!(aliquots_all %in% unique(fastq_info_df$aliquot))]
length(unique(fastq_info_df$file_name))
nrow(fastq_info_df)
## merge data 2
fastq_info_samplemap_df <- fastq_info_samplemap_df %>%
  select(-file_dir) %>%
  unique() %>%
  rename(file_name = `File Name`) %>%
  filter((data_type == "snRNA") & (file_name %in% fastq_info_df$file_name))
nrow(fastq_info_samplemap_df)
fastq_info_df = merge(x = fastq_info_df, y = fastq_info_samplemap_df, by = "file_name", all.x = T)
fastq_info_df$data_type = "snRNA"
fastq_info_df$`Flow Cell ID`[fastq_info_df$aliquot == "CPT0001220012"] = "HVLMGDSXX"
fastq_info_df$`Index Sequence`[fastq_info_df$aliquot == "CPT0001220012"] = "ATTACTTC-TGCGAACT-GCATTCGG-CAGCGGAA"
fastq_info_df$`Lane Number`[fastq_info_df$aliquot == "CPT0001220012"] = "4"
fastq_info_df$`Completion Date`[fastq_info_df$aliquot == "CPT0001220012"] = "10/12/19"
fastq_info_df$`Sample Name`[fastq_info_df$aliquot == "CPT0001220012"] = "TWCE-CPT000122-0012-lib1"
## merge data 3
snRNA_metadata_df <- merge(x = metadata_filtered_df, 
                           y = fastq_info_df, 
                           by.x = c("Aliquot.snRNA"), by.y = c("aliquot"), all.x = T)

# construct meta data -----------------------------------------------------------------
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
  mutate(Location =  paste0(str_split_fixed(file_path, "FASTQ\\/", 2)[,1], "FASTQ/")) %>%
  rename(`Sequencing lane` = `Lane Number`) %>%
  mutate(`Sample name` = ifelse(is.na(`Sample Name`), str_split_fixed(file_name, "\\-lib1", 2)[,1], `Sample Name`)) %>%
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
         Published, Status, Location, file_path, `Flow Cell ID`, `Index Sequence`, 
         `Sequencing lane`, `Sequence read length`, `Sample name`,
         `Number of expected cells`, `Sequencing date`, Protocol, md5sum)
snRNA_metadata_df$`Sequencing date`[snRNA_metadata_df$`Sequencing date` == "NA//"] <- "NA"
snRNA_metadata_df[is.na(snRNA_metadata_df)] <- "NA"

table(snRNA_metadata_df$snRNA_aliquot_id)
nrow(snRNA_metadata_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GDC.ccRCC.snRNA.FASTQ.", run_id, ".tsv")
write.table(x = snRNA_metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)


