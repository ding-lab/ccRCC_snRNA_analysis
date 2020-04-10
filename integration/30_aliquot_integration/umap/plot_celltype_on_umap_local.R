# Yige Wu @WashU Apr 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200406.v1/30_aliquot_integration.barcode2celltype.20200406.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)

# plot for cell group----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_group),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = cellgroup_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "cellgroup_on_umap.", run_id, ".png")
png(filename = file2write, width = 1100, height = 800, res = 150)
print(p)
dev.off()
# Tumor cells                  Immune                  Stroma Normal epithelial cells 
# "#7F3C8D"               "#11A579"               "#3969AC"               "#F2B701" 
# Unknown 
# "grey50" 

# plot for cell group----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type.shorter),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = celltype_shorter_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "celltype_shorter_on_umap.", run_id, ".png")
png(filename = file2write, width = 1600, height = 1000, res = 150)
print(p)
dev.off()
# Tumor cells                   Unknown                   B-cells                 Basophils 
# "#7F3C8D"                  "grey50"                 "#2E91E5"                 "#E15F99" 
# CD4/CD8 proliferating              CD4+ T-cells    CD8+ T-cells activated    CD8+ T-cells exhausted 
# "#1CA71C"                 "#FB0D0D"                 "#DA16FF"                 "#222A2A" 
# cDC               Macrophages           Macrophages M2b    Mixed myeloid/lymphoid 
# "#B68100"                 "#750D86"                 "#EB663B"                 "#511CFB" 
# NK cells                       pDC              Plasma cells                     Tregs 
# "#00A08B"                 "#FB00D1"                 "#FC0080"                 "#B2828D" 
# TRM         Endothelial cells               Fibroblasts            Myofibroblasts 
# "#6C7C32"                 "#778AAE"                 "#862A16"                 "#A777F1" 
# Distal convoluted tubule Intercalated cells type A Intercalated cells type B             Loop of Henle 
# "#87C55F"                 "#9EB9F3"                 "#FE88B1"                 "#C9DB74" 
# Podocytes           Proximal tubule 
# "#8BE0A4"                 "#B497E7" 

