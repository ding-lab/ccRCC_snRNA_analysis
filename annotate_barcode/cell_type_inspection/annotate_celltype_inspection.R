# Yige Wu @WashU Aug 2020

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
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input cell type inpection using mutation mapping
mutationmap_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/Inpection_Somatic_Mutation_Mapping.20200805.xlsx")

# map case ids and new aliquot ids ----------------------------------------
mutationmap_df$Id_Case <- mapvalues(x = mutationmap_df$Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
mutationmap_df$Id_Aliquot_WU <- mapvalues(x = mutationmap_df$Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
mutationmap_df <- mutationmap_df %>%
  arrange(Need_to_evaluate)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Inpection_Somatic_Mutation_Mapping.", run_id, ".tsv")
write.table(x = mutationmap_df, file = file2write, sep = "\t", quote = F, row.names = F)

