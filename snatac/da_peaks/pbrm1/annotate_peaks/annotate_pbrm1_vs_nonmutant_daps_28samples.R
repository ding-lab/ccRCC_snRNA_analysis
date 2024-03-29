# Yige Wu @WashU Jun 2021

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
## input all dap + caps
peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/UP_PBRM1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv")
peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/DOWN_PBRM1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv")
## input peak fold changes
peak2foldchange_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/PBRM1_vsNon_mutants.Filtered_peaks_byMinPct.0.1.20210712.tsv")
## input coaccessiblity results
coaccess_peak2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_coaccessible_peaks_28samples/20211011.v1/Coaccessible_Peaks.Annotated.20211011.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- rbind(peaks_up_df %>%
                         dplyr::rename(annotation = Type) %>%
                         dplyr::rename(Count_sig.snATAC = Count_up) %>%
                         dplyr::select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Up"),
                       peaks_down_df %>%
                         dplyr::rename(annotation = Type) %>%
                         dplyr::rename(Count_sig.snATAC = Count_down) %>%
                         dplyr::select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Down"))
peaks_anno_df$avg_log2FC <- mapvalues(x = peaks_anno_df$peak, from = peak2foldchange_df$peak, to = as.vector(peak2foldchange_df$avg_log2FC))
peaks_anno_df$avg_log2FC <- as.numeric(as.vector(peaks_anno_df$avg_log2FC))
peaks_anno_df <- peaks_anno_df %>%
  mutate(peak2gene_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1]) %>%
  select(-annotation)

# annotate with coaccessiblity results ------------------------------------
peaks_anno_wcoaccess_df <- merge(x = peaks_anno_df, 
                       y = coaccess_peak2genes_df %>%
                         select(Peak1, Peak2, peak2gene_type.2, genesymbol.2, coaccess) %>%
                         rename(peak.coaccess = Peak2) %>%
                         rename(peak2gene_type.coaccess = peak2gene_type.2) %>%
                         rename(genesymbol.coaccess = genesymbol.2) %>%
                         rename(coaccess_score = coaccess),
                       by.x = c("peak"), by.y = c("Peak1"))
peak2gene_enhancers_df <- peaks_anno_wcoaccess_df %>%
  filter(peak2gene_type != "Promoter" & peak2gene_type.coaccess == "Promoter")
peak2gene_enh_pro_df <- rbind(peak2gene_enhancers_df %>%
                                mutate(peak2gene_type = "Enhancer") %>%
                                mutate(Gene = genesymbol.coaccess),
                              peaks_anno_df %>%
                                filter(peak2gene_type == "Promoter") %>%
                                mutate(peak.coaccess = NA) %>%
                                mutate(peak2gene_type.coaccess = NA) %>%
                                mutate(genesymbol.coaccess = NA) %>%
                                mutate(coaccess_score = NA))

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_DAPs.Annotated.", run_id, ".tsv")
write.table(file = file2write, x = peaks_anno_df, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "PBRM1_DAP2Gene.EnhancerPromoter.", run_id, ".tsv")
write.table(file = file2write, x = peak2gene_enh_pro_df, quote = F, sep = "\t", row.names = F)


