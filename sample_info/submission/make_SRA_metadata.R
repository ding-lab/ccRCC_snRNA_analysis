# Yige Wu @ WashU 2021 Nov

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
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input clinical data
clinical_case_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data_snrna_samples/20211011.v1/snRNA_ccRCC_Clinicl_Table.20211011.v1.tsv")

# process -----------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, snRNA_available, snATAC_used, Sample_Type)
final_df <- rbind(metadata_filtered_df %>%
                    mutate(sample_name = paste0("snRNA.", Aliquot.snRNA.WU)) %>%
                    mutate(library_strategy = "RNA-seq") %>%
                    select(sample_name, Case, Sample_Type, library_strategy),
                  metadata_filtered_df %>%
                    filter(snATAC_used) %>%
                    mutate(sample_name = paste0("snATAC.", Aliquot.snRNA.WU)) %>%
                    mutate(library_strategy = "ATAC-seq") %>%
                    select(sample_name, Case, Sample_Type, library_strategy))
final_df <- merge(x = final_df, y = clinical_case_df, by = c("Case"), all.x = T)
final_df <- final_df %>%
  mutate(biosample_accession = "SAMN22866627") %>%
  mutate(library_ID = sample_name) %>%
  mutate(title = "Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma") %>%
  mutate(library_source = ifelse(library_strategy == "RNA-seq", "TRANSCRIPTOMIC SINGLE CELL", "OTHER")) %>%
  mutate(library_selection = "PCR") %>%
  mutate(library_layout = "Paired") %>%
  mutate(platform = "ILLUMINA") %>%
  mutate(instrument_model = "Illumina NovaSeq 6000") %>%
  mutate(design_description = paste0("Nuclei and barcoded beads were isolated in oil droplets via the 10x Genomics Chromium instrument. Single nuclei suspensions were counted and adjusted to a range of 500 to 1800 nuclei/µL using a hemocytometer. Reverse transcription was subsequently performed to incorporate cell and transcript specific barcodes.",
                                     ifelse(library_strategy == "RNA-seq", 
                                            " All snRNA-seq samples were run using the Chromium Next GEM Single Cell 3’ Library and Gel Bead Kit v3.1 (10x Genomics). Barcoded libraries were then pooled and sequenced on the Illumina NovaSeq 6000 system with specific flow cell type S4.",
                                            " For snATAC-seq, Chromium Next GEM Single Cell ATAC Library and Gel Bead Kit v1.1 prep (10x Genomics) were used for all samples. Barcoded libraries were then pooled and sequenced on the Illumina NovaSeq 6000 system with specific flow cell type S1")
                                     )) %>%
  mutate(filetype = "bam") %>%
  mutate(reference_fasta_file = "") %>%
  mutate(reference_assembly = "GRCh38") %>%
  mutate(filename = paste0(sample_name, ".bam")) %>%
  select(biosample_accession, library_ID, title, library_strategy, library_source, library_selection, library_layout, platform, instrument_model,
         design_description, filetype, reference_fasta_file, reference_assembly, filename)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "HumanTissue.SRA_Metadata.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


