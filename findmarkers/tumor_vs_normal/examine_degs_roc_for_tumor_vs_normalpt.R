# Yige Wu @WashU Apr 2020

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
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## input DEGs
degfile_names <- list.files(path = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/", recursive = T)
degfile_names
degfile_names <- degfile_names[grepl(pattern = "findmarker_roc", x = degfile_names)]
degfile_names

# read deg files ----------------------------------------------------------
deg_df <- NULL
for (degfile_name in degfile_names) {
  path_degfile <- paste0("./Resources/Analysis_Results/findmarkers/tumor_vs_normal/", degfile_name)
  deg_tmp <- fread(input = path_degfile)
  deg_tmp <- deg_tmp %>%
    mutate(aliquot_tumor = str_split_fixed(string = degfile_name, pattern = "_", n = 4)[,3]) %>%
    rename(gene_symbol = row_name)
  deg_df <- rbind(deg_df, deg_tmp)
}
## classify genes in to up or down regulated in tumor vs normal proximal tubule
deg_df <- deg_df %>%
  mutate(direction_tumor_vs_normal = ifelse(avg_diff > 0, "tumor+", "tumor-"))
## make summary of genes
deg_summary_df <- deg_df %>%
  group_by(gene_symbol, direction_tumor_vs_normal) %>%
  summarise(Freq = n(), mean_avg_diff = mean(avg_diff), mean_pct.1 = mean(pct.1), mean_pct.2 = mean(pct.2), mean_power = mean(power)) %>%
  arrange(-Freq, -mean_power)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findmarkers_roc.", "tumor_vs_normalpt.", run_id, ".tsv")
write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
file2write <- paste0(dir_out, "summary_findmarkers_roc.", "tumor_vs_normalpt.", run_id, ".tsv")
write.table(x = deg_summary_df, file = file2write, sep = "\t", quote = F, row.names = F)

# examine known pathway members -------------------------------------------
## input the hif pathway members
hiftargets_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200302.v1/HIF_Target_Genes.20200302.v1.tsv", data.table = F)
deg_summary_hif_df <- deg_summary_df %>%
  filter(gene_symbol %in% hiftargets_df$target_genesymbol)
## input the hif pathway members
deg_summary_epigeneticregulators_df <- deg_summary_df %>%
  filter(gene_symbol %in% c(pbaf_genes, "SETD2", "BAP1", "KDM5C"))
## input the hif pathway members
load(file = "./Resources/Gene_Lists/2015-08-01_Gene_Set.RData")
deg_summary_pi3kmtor_df <- deg_summary_df %>%
  filter(gene_symbol %in% REACT[["165159	mTOR signalling"]])

