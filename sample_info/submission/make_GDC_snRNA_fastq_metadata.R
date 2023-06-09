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
# fastq_summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.tsv")
fastq_summary_list <- readRDS("./Resources/Analysis_Results/sample_info/submission/extract_fastq_summary/20221212.v1/ccRCC.snRNA.snATAC.fastq.summary.20221212.v1.RSD")

# preprocess fastq summary info -------------------------------------------
fastqinfo_summary_df <- NULL
for (idx in names(fastq_summary_list)) {
  fastqinfo_tmp_df <- fastq_summary_list[[idx]]
  if ("snRNA" %in% fastqinfo_tmp_df$data_type) { ## we are only extracting the snRNA FASTQ files
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
  mutate(file_path = paste0("/diskmnt/primary/ccRCC_snRNA/", Aliquot.snRNA, "/BAM/possorted_genome_bam.bam")) %>%
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

# combine -----------------------------------------------------------------
final_df <- rbind(snRNA_metadata_df, snATAC_metadata_df)
final_df$`Sequencing date`[final_df$`Sequencing date` == "NA//"] <- NA

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GDC.ccRCC.snRNA.FASTQ.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


