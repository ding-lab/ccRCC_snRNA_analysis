# Yige Wu @WashU Feb 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution is too confusing

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
### set output subdirectories
dir_out1 <- paste0(dir_out, "All_CNVs/")
dir.create(dir_out1)
### set output subdirectories
dir_out2 <- paste0(dir_out, "Expected_CNVs/")
dir.create(dir_out2)
### set output subdirectories
dir_out3 <- paste0(dir_out, "Expected_CNVs_Higher_Than_NAT/")
dir.create(dir_out3)


# input dependencies ------------------------------------------------------
## load meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_state_count_aliquots <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_malignant_nephron_epithelium_per_cluster/20200225.v1/Fraction_of_Malignant_Nephron_Epithelium_with_CNA_by_Gene.Per_Tumor_Subcluster.20200225.v1.tsv", data.table = F)
table(cnv_state_count_aliquots$tumor_subcluster)
## load CNV fraction in NAT sample
cnv_state_count_nat <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_nephron_epithelium_per_sample/20200226.v1/Fraction_of_Nephron_Epithelium_with_CNA_by_Gene.20200226.v1.tsv", data.table = F)
cnv_state_count_nat <- cnv_state_count_nat %>%
  filter(aliquot %in% meta_tab$Aliquot.snRNA[meta_tab$Sample_Type == "Normal"])


for (aliquot_tmp in unique(cnv_state_count_aliquots$aliquot)) {
  # make data frame for plotting --------------------------------------------
  plot_data_df <- cnv_state_count_aliquots %>%
    filter(aliquot == aliquot_tmp)
  ## filter out copy neutral ones
  plot_data_df <- plot_data_df %>%
    filter(cna_state != 1)
  ## add text for cnv state to show
  plot_data_df$sn_cna_cat <- map_infercnv_state2category(copy_state = plot_data_df$cna_state)
  ## order chromosome region
  plot_data_df$chr_region <- factor(plot_data_df$chr_region, levels = c("3p", "5q", "14q",
                                                                        "9p21", 
                                                                        "10q23",
                                                                        "1p31",
                                                                        "6q24",
                                                                        "3p12",
                                                                        "9p23",
                                                                        "14q24",
                                                                        "3p26",
                                                                        "1q32",
                                                                        "8q24",
                                                                        "9q24"))
  ## order the fraction
  plot_data_df <- plot_data_df %>%
    arrange(desc(Fraction))
  ## make cluster id into factor
  plot_data_df$tumor_subcluster <- as.factor(plot_data_df$tumor_subcluster)
  ## get case id
  case_tmp <- meta_tab$Case[meta_tab$Aliquot.snRNA == aliquot_tmp]
  
  # plot expected + unexpected CNVs--------------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, 
                      mapping = aes(x = gene_symbol, y = tumor_subcluster, size = Fraction, fill = sn_cna_cat, color = expected_cna), 
                      shape = 21, alpha = 0.8)
  p <- p + scale_fill_manual(values = copy_number_colors)
  p <- p + scale_color_manual(values = c("TRUE" = "#000000", "FALSE" = "#FFFFFF"))
  p <- p + facet_grid(.~chr_region + gene_cna_type, scales = "free", space = "free", shrink = T)
  p <- p + ggtitle(label = paste0("CNV Distribution in ", aliquot_tmp, " Tumor Subclusters"))
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12, face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(panel.spacing = unit(0, "lines"))
  p <- p + theme(strip.text.y = element_text(angle = 0))
  p
  png(file = paste0(dir_out1, case_tmp, ".", aliquot_tmp, ".", "Fraction_of_Malignant_Nephron_Epithelium_with_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".", run_id, ".png"), 
      width = 1000, height = 600, res = 150)
  print(p)
  dev.off()
  
  # plot expected CNVs--------------------------------------------------------------------
  plot_data_df <- plot_data_df %>%
    filter(expected_cna)
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, 
                      mapping = aes(x = gene_symbol, y = tumor_subcluster, size = Fraction, fill = sn_cna_cat, color = expected_cna), 
                      shape = 21, alpha = 0.8)
  p <- p + scale_fill_manual(values = copy_number_colors)
  p <- p + scale_color_manual(values = c("TRUE" = "#000000", "FALSE" = "#FFFFFF"))
  p <- p + facet_grid(.~chr_region + gene_cna_type, scales = "free", space = "free", shrink = T)
  p <- p + ggtitle(label = paste0("CNV Distribution in ", aliquot_tmp, " Tumor Subclusters"))
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12, face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(panel.spacing = unit(0, "lines"))
  p <- p + theme(strip.text.y = element_text(angle = 0))
  p
  png(file = paste0(dir_out2, case_tmp, ".", aliquot_tmp, ".", "Fraction_of_Malignant_Nephron_Epithelium_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".", run_id, ".png"), 
      width = 1000, height = 600, res = 150)
  print(p)
  dev.off()
  
  
  # plot expected CNVs but higher than the fraction in the NAT sample--------------------------------------------------------------------
  plot_data_df <- plot_data_df %>%
    filter(expected_cna)
  ## add in the fraction for CNV in NAT (same gene, same cnv state)
  plot_data_df <- base::merge(plot_data_df, cnv_state_count_nat[, c("gene_symbol", "cna_state", "Fraction")], 
                              by = c("gene_symbol", "cna_state"), suffix = c("", ".nat"), all.x = T)
  ## make NA value 0
  tmp <- plot_data_df$Fraction.nat
  tmp[is.na(tmp)] <- 0
  plot_data_df$Fraction.nat <- tmp
  ## filter
  plot_data_df <- plot_data_df %>%
    filter(Fraction > Fraction.nat)
  ## order the fraction
  plot_data_df <- plot_data_df %>%
    arrange(desc(Fraction))
  
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, 
                      mapping = aes(x = gene_symbol, y = tumor_subcluster, size = Fraction, fill = sn_cna_cat, color = expected_cna), 
                      shape = 21, alpha = 0.8)
  p <- p + scale_fill_manual(values = copy_number_colors)
  p <- p + scale_color_manual(values = c("TRUE" = "#000000", "FALSE" = "#FFFFFF"))
  p <- p + facet_grid(.~chr_region + gene_cna_type, scales = "free", space = "free", shrink = T)
  p <- p + ggtitle(label = paste0("CNV Distribution in ", aliquot_tmp, " Tumor Subclusters"))
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12, face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(panel.spacing = unit(0, "lines"))
  p <- p + theme(strip.text.y = element_text(angle = 0))
  p
  png(file = paste0(dir_out3, case_tmp, ".", aliquot_tmp, ".", "Fraction_of_Malignant_Nephron_Epithelium_with_Expected_CNA_by_Gene", 
                    ".", "Higher_than_NAT", ".", "By_Tumor_Subcluster", ".", run_id, ".png"), 
      width = 1000, height = 600, res = 150)
  print(p)
  dev.off()
}
