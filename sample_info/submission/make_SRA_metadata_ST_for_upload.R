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

# process -----------------------------------------------------------------
final_df <- data.frame(library_ID = c("HT282N1_R1", "HT282N1_R2", "HT293N1_R1", "HT293N1_R2"),
                       biosample_accession = c("SAMN29632060", "SAMN29632060", "SAMN29632061", "SAMN29632061"),
                       filename = c("TWCE-HT282N1-S1H3Fs4U1Bp1-HT282N1-S1H3Fs4U1Bp1_2_S1_L004_R1_001.fastq.gz",
                                    "TWCE-HT282N1-S1H3Fs4U1Bp1-HT282N1-S1H3Fs4U1Bp1_2_S1_L004_R2_001.fastq.gz",
                                    "TWCE-HT293N1-S1H3Fs1U1Bp1-HT293N1-S1H3Fs1U1Bp1_2_S4_L004_R1_001.fastq.gz",
                                    "TWCE-HT293N1-S1H3Fs1U1Bp1-HT293N1-S1H3Fs1U1Bp1_2_S4_L004_R2_001.fastq.gz"))
final_df <- final_df %>%
  mutate(title = "Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma") %>%
  mutate(library_strategy = "RNA-seq") %>%
  mutate(library_source = "TRANSCRIPTOMIC SINGLE CELL") %>%
  mutate(library_selection = "PCR") %>%
  mutate(library_layout = "Paired") %>%
  mutate(platform = "ILLUMINA") %>%
  mutate(instrument_model = "Illumina NovaSeq 6000") %>%
  mutate(design_description = "The RNA quality of FFPE tissue blocks was evaluated by calculating DV200 of RNA extracted from FFPE tissue sections following the Qiagen RNeasy FFPE Kit protocol. After the Tissue Adhesion Test, 5 μm sections were placed on the Visium Spatial Gene Expression Slide following Visium Spatial Protocols-Tissue Preparation Guide (10X Genomics, CG000408 Rev A). After overnight drying, slides were incubated at 60°C for 2 h. Deparaffinization was then performed following  Visium Spatial for FFPE – Deparaffinization, H&E Staining, Imaging & Decrosslinking Protocol (10X Genomics, CG000409 Rev A). Sections were stained with hematoxylin and eosin and imaged at 20x magnification using the brightfield imaging setting on a Leica DMi8 microscope. After that, decrosslinking was performed immediately for H&E stained sections. Next, human whole transcriptome probe panels were then added to the tissue. After these probe pairs hybridized to their target genes and ligated to one another, the ligation products were released following RNase treatment and permeabilization. The ligated probes were then hybridized to the spatially barcoded oligonucleotides on the Capture Area. Spatial Transcriptomics libraries were generated from the probes and sequenced on the S4 flow cell of the Illumina NovaSeq 6000 system.") %>%
  mutate(filetype = "fastq") %>%
  mutate(path_katmai = paste0("/diskmnt/primary/Spatial_Transcriptomics/data/globus/2864517/", filename)) %>%
  arrange(library_strategy) %>%
  select(biosample_accession, library_ID, title, library_strategy, library_source, library_selection, library_layout, platform, instrument_model,
         design_description, filetype, filename, path_katmai)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ST.HumanTissue.SRA_Metadata.forupload.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


