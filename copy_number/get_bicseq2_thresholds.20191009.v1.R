# Yige Wu @ WashU Oct 2019
# compare BICSEQ2 with TCGA result from cbioportal

# https://www.cbioportal.org/faq#what-do-amplification-gain-deep-deletion-shallow-deletion-and--2--1-0-1-and-2-mean-in-the-copy-number-data
# What do “Amplification”, “Gain”, “Deep Deletion”, “Shallow Deletion” and "-2", "-1", "0", "1", and "2" mean in the copy-number data? 
#   These levels are derived from copy-number analysis algorithms like GISTIC or RAE, and indicate the copy-number level per gene:
#   
#   -2 or Deep Deletion indicates a deep loss, possibly a homozygous deletion
# -1 or Shallow Deletion indicates a shallow loss, possibley a heterozygous deletion
# 0 is diploid
# 1 or Gain indicates a low-level gain (a few additional copies, often broad)
# 2 or Amplification indicate a high-level amplification (more copies, often focal)
# Note that these calls are putative. We consider the deep deletions and amplifications as biologically relevant for individual genes by default. Note that these calls are usually not manually reviewed, and due to differences in purity and ploidy between samples, there may be false positives and false negatives.
# 

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# functions ---------------------------------------------------------------
plot_tcga_cptac_pct <- function(plot_df, this_cnv_type, id, this_cnv_threshold, labels = FALSE, dir_out){
  p <- plot_df %>% dplyr::filter(cnv_type == this_cnv_type) %>%
    ggplot(aes(x = TCGA_pct, y = CPTAC_pct, label = Gene)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_smooth(method = "lm") +
    geom_point() +
    xlim(0,100) + ylim(0, 100) +
    labs(x = "TCGA Percentage with CNV Event",
         y = "CPTAC Percentage with CNV Event",
         title = str_c("TCGA/CPTAC Gene Level CNV Comparison\n(",
                       this_cnv_type, ": log2(copy_ratio) = ", signif(x = this_cnv_threshold, digits = 3), ")")) +
    theme_bw(base_size = 20)
  if (labels) {
    p <- p + geom_label_repel()
    filename = str_c(dir_out, "TCGA-CPTAC_CNV_comparision.", id, ".", this_cnv_type, "_cutoff", signif(x = this_cnv_threshold, digits = 3), ".with_labels.png")
  } else {
    filename = str_c(dir_out, "TCGA-CPTAC_CNV_comparision.", id, ".", this_cnv_type, "_cutoff", signif(x = this_cnv_threshold, digits = 3), "without_labels.png")
  }
  png(filename, width = 600, height = 480)
  print(p)
  dev.off()
}

tcga_pct <- function(tcga_tbl, gene_name, cnv_type){
  gene_values <- tcga_tbl %>% dplyr::filter(Hugo_Symbol == gene_name) %>% dplyr::select(-Hugo_Symbol) %>% dplyr::select(-Entrez_Gene_Id)
  if (cnv_type == "Amplification"){
    return_pct <- mean(gene_values >= 1, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values <= -1, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values == 0, na.rm = T)*100
  }
  return(return_pct)
}

## get function for calculating TCGA
cptac_pct <- function(cptac_tbl, gene_name, cnv_type,
                      amp_threshold, del_threshold){
  gene_values <- cptac_tbl %>% dplyr::filter(gene == gene_name) %>% dplyr::select(-gene)
  if (cnv_type == "Amplification") {
    return_pct <- mean(gene_values > amp_threshold, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values < del_threshold, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values >= del_threshold & gene_values <= amp_threshold, na.rm = T)*100
  }
  return(return_pct)
}


tcga_pct_deep <- function(tcga_tbl, gene_name, cnv_type){
  gene_values <- tcga_tbl %>% dplyr::filter(Hugo_Symbol == gene_name) %>% dplyr::select(-Hugo_Symbol) %>% dplyr::select(-Entrez_Gene_Id)
  if (cnv_type == "Amplification"){
    return_pct <- mean(gene_values >= 2, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values <= -2, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values == 0, na.rm = T)*100
  }
  return(return_pct)
}


# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input ----------------------------------------------------------------
cancer <- "ccRCC"
cptac_tbl <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/copy_number/parse_bicseq2_cnv/20191009.v1/CPTAC3_ccRCC_Discovery_Set_BICSEQ2.20191009.v1.tsv", data.table = F)
tcga_tbl <- fread(file = "./Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pub/data_CNA.txt", data.table = F)

# get percentage of shallow amp & del -------------------------------------
## get the percentage of cBioPortal CNV (top)
top_amp_genes <- tcga_tbl$Hugo_Symbol[rowSums(tcga_tbl >=1) >= sort(rowSums(tcga_tbl >= 1), decreasing = T)[50]]
top_del_genes <- tcga_tbl$Hugo_Symbol[rowSums(tcga_tbl <= -1) >= sort(rowSums(tcga_tbl <= -1), decreasing = T)[50]]

cptac_gene_list <- cptac_tbl %>% pull(gene)
tcga_gene_list <- c(top_amp_genes, top_del_genes)


shared_gene_list <- tcga_gene_list[tcga_gene_list %in% cptac_gene_list]
n_shared_genes <- length(shared_gene_list)

# tcga_tbl <- as.tibble(tcga_tbl)
# cptac_tbl <- as.tibble(cptac_tbl)

for (del_threshold in c(-0.2, -0.1, -0.05)) {
  plot_df <- tibble(Gene = rep(shared_gene_list, 3),
                    cnv_type = c(rep("Amplification", n_shared_genes),
                                 rep("Deletion", n_shared_genes),
                                 rep("Neutral", n_shared_genes)),
                    TCGA_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Amplification")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Deletion")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Neutral"))),
                    CPTAC_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Amplification", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Deletion", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Neutral", amp_threshold, del_threshold))))
  plot_tcga_cptac_pct(plot_df, "Deletion", paste0("shallow", ".", cancer), this_cnv_threshold = del_threshold, labels = T, dir_out = dir_out)
  plot_tcga_cptac_pct(plot_df, "Deletion", paste0("shallow", ".", cancer), this_cnv_threshold = del_threshold, labels = F, dir_out = dir_out)
  
}
print(tcga_pct(tcga_tbl, "TP53", "Deletion"))
print(cptac_pct(cptac_tbl, "TP53", "Deletion", log2(1.1), log2(0.9)))

for (amp_threshold in c(0.2, 0.1, 0.05)) {
  # for (amp_threshold in c(log2(1.1))) {
  plot_df <- tibble(Gene = rep(shared_gene_list, 3),
                    cnv_type = c(rep("Amplification", n_shared_genes),
                                 rep("Deletion", n_shared_genes),
                                 rep("Neutral", n_shared_genes)),
                    TCGA_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Amplification")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Deletion")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct(tcga_tbl, x, "Neutral"))),
                    CPTAC_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Amplification", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Deletion", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Neutral", amp_threshold, del_threshold))))
  plot_tcga_cptac_pct(plot_df, "Amplification", paste0("shallow", ".", cancer), this_cnv_threshold = amp_threshold, labels = T, dir_out = dir_out)
  plot_tcga_cptac_pct(plot_df, "Amplification", paste0("shallow", ".", cancer), this_cnv_threshold = amp_threshold, labels = F, dir_out = dir_out)
  
}

next()
# get percentage of deep amp & del -------------------------------------
## get the percentage of cBioPortal CNV (top)
top_amp_genes <- tcga_tbl$Hugo_Symbol[rowSums(tcga_tbl >= 2) >= sort(rowSums(tcga_tbl >= 2), decreasing = T)[100]]
top_del_genes <- tcga_tbl$Hugo_Symbol[rowSums(tcga_tbl <= -2) >= sort(rowSums(tcga_tbl <= -2), decreasing = T)[100]]

cptac_gene_list <- cptac_tbl %>% pull(gene)
tcga_gene_list <- c(top_amp_genes, top_del_genes)


shared_gene_list <- tcga_gene_list[tcga_gene_list %in% cptac_gene_list]
n_shared_genes <- length(shared_gene_list)

for (amp_threshold in c(0.3, 0.4, 0.5, 0.6)) {
  plot_df <- tibble(Gene = rep(shared_gene_list, 3),
                    cnv_type = c(rep("Amplification", n_shared_genes),
                                 rep("Deletion", n_shared_genes),
                                 rep("Neutral", n_shared_genes)),
                    TCGA_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Amplification")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Deletion")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Neutral"))),
                    CPTAC_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Amplification", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Deletion", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Neutral", amp_threshold, del_threshold))))
  plot_tcga_cptac_pct(plot_df, "Amplification", "deep", this_cnv_threshold = amp_threshold, dir_out = dir_out)
}

for (del_threshold in c(-0.3, -0.4, -0.5, -0.6)) {
  plot_df <- tibble(Gene = rep(shared_gene_list, 3),
                    cnv_type = c(rep("Amplification", n_shared_genes),
                                 rep("Deletion", n_shared_genes),
                                 rep("Neutral", n_shared_genes)),
                    TCGA_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Amplification")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Deletion")),
                                 apply(as.matrix(shared_gene_list), 1, function(x) tcga_pct_deep(tcga_tbl, x, "Neutral"))),
                    CPTAC_pct = c(apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Amplification", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Deletion", amp_threshold, del_threshold)),
                                  apply(as.matrix(shared_gene_list), 1, function(x) cptac_pct(cptac_tbl, x, "Neutral", amp_threshold, del_threshold))))
  plot_tcga_cptac_pct(plot_df, "Deletion", "deep", this_cnv_threshold = del_threshold, dir_out = dir_out)
}
