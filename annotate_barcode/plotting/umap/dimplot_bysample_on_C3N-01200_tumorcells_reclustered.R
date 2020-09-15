# Yige Wu @WashU Sep 2020

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
## input the integrated data
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-01200_tumorlikecells_and_recluster/20200910.v1/TumorLikeCells.Reclustered.20200910.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# Dimplot -----------------------------------------------------------------
## map aliquot id
srat@meta.data$id_aliquot <- mapvalues(x = srat@meta.data$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## make distinguishable colors
colors_cluster <- Polychrome::dark.colors(n = length(unique(srat@meta.data$id_aliquot)))
names(colors_cluster) <- unique(srat@meta.data$id_aliquot)
p <- DimPlot(object = srat, group.by = "id_aliquot")
p <- p + scale_color_manual(values = colors_cluster)
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Sample.", "C3N-01200", ".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()

