# Yige Wu @WashU March 2020
## unite all the 10xmapping result

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
## set 10XMapping output directory
dir_10xmapping <- paste0("./Resources/snRNA_Processed_Data/10Xmapping/outputs/run_samples.snATAC.katmai/")

# input the 10xmapping result from each sample ----------------------------
mutation_map_sup_tab <- NULL
for (snRNA_aliquot_id_tmp in list.files(dir_10xmapping)) {
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab <- fread(input = paste0(dir_10xmapping, 
                                           snRNA_aliquot_id_tmp, "/", 
                                           snRNA_aliquot_id_tmp, "_mapping_heatmap_0.txt"), data.table = F)
  mutation_map_tab <- mutation_map_tab %>%
    mutate(gene_symbol = str_split_fixed(string = Mutation, pattern = "-", n = 3)[,1])
  mutation_map_tab.m <- melt(mutation_map_tab, id.vars = c("gene_symbol", "Mutation"))
  mutation_map_tab.m <- mutation_map_tab.m %>%
    mutate(allele_type = str_split_fixed(string = Mutation, pattern = "-", n = 3)[,3])
  
  mutation_map_tab.filtered <- mutation_map_tab.m %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable) %>%
    rename(mutation = Mutation) %>%
    mutate(aliquot = snRNA_aliquot_id_tmp)
  
  mutation_map_sup_tab <- rbind(mutation_map_tab.filtered, mutation_map_sup_tab)
}

# write table -------------------------------------------------------------
write.table(x = mutation_map_sup_tab, file = paste0(dir_out, "10XMapping.snATAC.mapping_heatmap.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

