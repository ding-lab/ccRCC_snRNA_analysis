# Yige Wu @WashU Aug 2021

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
deg_ccRCC_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")
deg_sex_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/PT_sex_DEGs.xlsx")

# overlap -----------------------------------------------------------------
deg_sex_df$genesymbol_human <- toupper(x = deg_sex_df$`PT sex DEG list`)
deg_filtered_df <- deg_ccRCC_df %>%
  filter(genesymbol_deg %in% deg_sex_df$genesymbol_human)
table(deg_filtered_df$Tumor_vs_PT)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT.snRNA_DEGs.OverlapPTSexDEGs.", run_id, ".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
