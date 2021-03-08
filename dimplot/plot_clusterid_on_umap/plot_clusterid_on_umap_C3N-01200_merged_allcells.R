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
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200828.v1/31Aliquot.Barcode2CellType.20200828.v1.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample//20200406.v1/individual_cluster_meta_data.20200406.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# input seurat object, and get umap info -----------------------------------------------------
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")
umap_df <- FetchData(object = srat, vars = c("orig.ident", "UMAP_1", "UMAP_2", "seurat_clusters"))
umap_df$barcode_integrated <- rownames(umap_df)


# plot for cell group----------------------------------------------------------
aliquot_show <- "C3N-01200 Merged"
plotdata_df <- umap_df
colors_cluster <- Polychrome::alphabet.colors(n = 18)
names(colors_cluster) <- 0:17
p <- ggplot()
p <- p + geom_point(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = colors_cluster)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + ggtitle(label = paste0(aliquot_show, " Clusters"))
p
file2write <- paste0(dir_out, aliquot_show, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()
