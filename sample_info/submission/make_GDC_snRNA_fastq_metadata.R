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

# preprocess meta data --------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, Sample, snRNA_available, snATAC_used, Sample_Type, Aliquot.snRNA)
aliquots_all = unique(metadata_filtered_df$Aliquot.snRNA)

# preprocess fastq summary info -------------------------------------------
fastq_info_samplemap_df <- NULL
for (idx in names(fastq_summary_list)) {
  fastqinfo_tmp_df <- fastq_summary_list[[idx]]
  if ("snRNA" %in% fastqinfo_tmp_df$data_type) { ## we are only extracting the snATAC FASTQ files
    colname_samplename <- colnames(fastqinfo_tmp_df)[grepl(pattern = "Sample", x = colnames(fastqinfo_tmp_df))]
    colname_lane <- colnames(fastqinfo_tmp_df)[grepl(pattern = "Lane", x = colnames(fastqinfo_tmp_df))]
    if ("File Name" %in% colnames(fastqinfo_tmp_df)) {
      fastqinfo_keep_tmp_df <- fastqinfo_tmp_df[, c("File Name", "Flow Cell ID", "Index Sequence", "Completion Date", colname_samplename, colname_lane)]
    } else {
      ## in this case, there are two FASTQ files listed in each row
      colnames2process <- colnames(fastqinfo_tmp_df)[which(grepl(x = fastqinfo_tmp_df[1,], pattern = "fastq"))]
      fastqinfo_tmp_df1 <- fastqinfo_tmp_df[, c("Flow Cell ID", "Index Sequence", "Completion Date", colname_samplename, colname_lane, colnames2process)]
      fastqinfo_tmp_df2 <- melt.data.table(data = data.table(fastqinfo_tmp_df1), measure.vars = colnames2process)
      fastqinfo_tmp_df2 <- as.data.frame(fastqinfo_tmp_df2)
      fastqinfo_keep_tmp_df <- fastqinfo_tmp_df2[, c( "value", "Flow Cell ID", "Index Sequence", "Completion Date", colname_samplename, colname_lane)]
    }
    colnames(fastqinfo_keep_tmp_df) <- c("File Name", "Flow Cell ID", "Index Sequence", "Completion Date", "Sample Name", "Lane Number")
    fastqinfo_keep_tmp_df$`Completion Date` = as.character(fastqinfo_keep_tmp_df$`Completion Date`)
    fastq_info_samplemap_df <- rbind(fastq_info_samplemap_df, fastqinfo_keep_tmp_df)
  }
}

## process MD5
md5_df = md5_df %>%
  rename(file_name.md5 = file_name) %>%
  mutate(file_name = gsub(x = file_name.md5, pattern = "\\.md5", replacement = "")) %>%
  rename(file_path.md5 = file_path) %>%
  mutate(file_path = gsub(x = file_path.md5, pattern = "\\.md5", replacement = "")) %>%
  mutate(file_path = gsub(x = file_path, pattern = "\\/MD5\\/", replacement = "/FASTQ/"))
## merge data
fastq_info_df = merge(x = fastq_paths_df, y = md5_df, by = "file_path")
fastq_info_df = fastq_info_df %>%
  filter(grepl(x = file_path, pattern = "ccRCC_snRNA")) %>%
  mutate(aliquot = str_split_fixed(file_path, "\\/", 6)[,5]) %>%
  filter(aliquot %in% metadata_filtered_df$Aliquot.snRNA)
length(unique(fastq_info_df$aliquot))
aliquots_all[!(aliquots_all %in% unique(fastq_info_df$aliquot))]
fastq_info_df = merge(x = fastq_info_df, y = fastq_info_samplemap_df, by.x = "file_name", by.y = "File Name", all.x = T)

# process snRNA data-----------------------------------------------------------------
snRNA_metadata_df <- merge(x = metadata_filtered_df, 
                            y = fastq_info_df, 
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
  mutate(Location =  paste0(str_split_fixed(file_path, "FASTQ\\/", 2)[,1], "FASTQ/")) %>%
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
         Published, Status, Location, file_path, `Flow Cell ID`, `Index Sequence`, 
         `Sequencing lane`, `Sequence read length`, `Sample name`,
         `Number of expected cells`, `Sequencing date`, Protocol)

# combine -----------------------------------------------------------------
final_df <- rbind(snRNA_metadata_df, snATAC_metadata_df)
final_df$`Sequencing date`[final_df$`Sequencing date` == "NA//"] <- NA

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GDC.ccRCC.snRNA.FASTQ.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


