# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20201207.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20201207.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_top_geneset_scores/20210419.v1/MSigDB.Hallmark.tsv")

# set parameters ----------------------------------------------------------
## specify the CNV gene and type to plot
cnvs_plot_df <- knowncnvgenes_df
## specify the pathway scores to plot
genesetnames_plot <- c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_ALLOGRAFT_REJECTION",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING",
                       "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_IL2_STAT5_SIGNALING",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_ESTROGEN_RESPONSE_EARLY")

# make data frame for plotting --------------------------------------------
## preprocess
scores_df <- scores_df %>%
  mutate(tumor_subcluster = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
i <- 31
for (genesetname in genesetnames_plot[2:length(genesetnames_plot)]) {
  ## extract pathway score
  geneset_colname <- paste0(gsub(x = genesetname, pattern = "HALLMARK_", replacement = ""), "_Score")
  scores_tmp_df <- scores_df[, c("tumor_subcluster", geneset_colname)] 
  colnames(scores_tmp_df) <- c("tumor_subcluster", "geneset_score")
  dir_out1 <- paste0(dir_out, genesetname, "/"); dir.create(dir_out1)
  for (i in 1:nrow(cnvs_plot_df)) {
    ## extract CNV data
    gene_cnv <- cnvs_plot_df$Gene_Symbol[i]
    cnv_type <- cnvs_plot_df$CNV_Type[i]
    cnv_tmp_df <- cnv_3state_count_aliquots %>%
      filter(gene_symbol == gene_cnv) %>%
      filter(cna_3state == cnv_type)
    if (nrow(cnv_tmp_df) == 0) {
      next()
    }
    ## make data frame for the NA data
    cluster_nonna_df <- cnv_3state_count_aliquots %>%
      filter(gene_symbol %in% gene_cnv) %>%
      select(tumor_subcluster) %>%
      unique() %>%
      mutate(Data_detected = T)
    ## merge
    plotdata_df <- merge(x = cnv_tmp_df, y = scores_tmp_df, by = c("tumor_subcluster"), all.y = T)
    plotdata_df <- merge(x = plotdata_df, y = cluster_nonna_df, by = c("tumor_subcluster"), all.x = T)
    ## filter
    plotdata_df <- plotdata_df %>%
      filter(!is.na(Data_detected)) %>%
      mutate(Fraction = ifelse(is.na(Fraction), 0, Fraction)) %>%
      mutate(sample_id = str_split_fixed(string = tumor_subcluster, pattern = "_", n = 2)[,1])
    
    ## plot
    p <- ggplot()
    p <- p + geom_line(data = plotdata_df, mapping = aes(x = Fraction, y = geneset_score, group = sample_id))
    p <- p + geom_point(data = plotdata_df, mapping = aes(x = Fraction, y = geneset_score), alpha = 0.8)
    p <- p + theme_classic()
    p <- p + facet_wrap(facets = c("sample_id"), nrow = 3)
    p <- p + ylab(label = geneset_colname) + xlab(paste0("% ", gene_cnv, " ", cnv_type, " per tumor cluster"))
    file2write <- paste0(dir_out1, gene_cnv, "_", cnv_type, ".", genesetname, ".png")
    png(file2write, width = 1500, height = 800, res = 150)
    print(p)
    dev.off()
  }
}




