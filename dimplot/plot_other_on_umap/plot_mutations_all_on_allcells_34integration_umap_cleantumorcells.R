# Yige Wu @WashU Jan 2022
## plot mutation mapping type on integration UMAP
## clean out tumor cells in non-tumor clusters

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20220126.v1/10XMapping.20220126.v1.tsv", data.table = F)
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)

# make plot data----------------------------------------------------------
integrated_umap_df <- integrated_umap_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(id_barcode = paste0(orig.ident, "_", barcode_individual))
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(id_barcode = paste0(aliquot, "_", barcode)) %>%
  mutate(id_mutation = gsub(x = mutation, pattern = "\\-Var||\\-Ref", replacement = "")) %>%
  mutate(id_aliquot_mut = paste0(aliquot, "_", id_mutation))
barcodes_mut_df <- snRNA_mutation_df %>%
  filter(allele_type == "Var")
barcodes_ref_df <- snRNA_mutation_df %>%
  filter(allele_type == "Ref") %>%
  filter(id_aliquot_mut %in% barcodes_mut_df$id_aliquot_mut)
plot_data_df <- integrated_umap_df %>%
  mutate(mutation_mapped_type = ifelse(id_barcode %in% barcodes_mut_df$id_barcode,
                                       "Var",
                                       ifelse(id_barcode %in% barcodes_ref_df$id_barcode, "Ref", "NA")))
plot_data_df <- merge(plot_data_df,
                      barcode2celltype_df %>%
                        mutate(Cell_group = Cell_group_w_epithelialcelltypes) %>%
                        select(orig.ident, individual_barcode, Cell_group),
                      by.x = c("orig.ident", "barcode_individual"),
                      by.y = c("orig.ident", "individual_barcode"), all.x = T
)
# table(plot_data_df$mutation_mapped_type)
plot_data_df$ident <- as.character(plot_data_df$ident)
## consider the tumor cells in non-tumor clusters unknown and remove
clusterids_nontumor <- c("2", "29", "32")
plot_data_df <- plot_data_df %>%
  mutate(is_unknown = (ident %in% clusterids_nontumor) & (Cell_group == "Tumor cells")) %>%
  filter(!is_unknown)
plot_data_df %>%
  filter(ident %in% clusterids_nontumor) %>%
  filter(mutation_mapped_type == "Var") %>%
  select(Cell_group) %>%
  table()

### order the data frame so that cells mapped with variant allele will come on top
plot_data_df <- rbind(plot_data_df[plot_data_df$mutation_mapped_type == "NA",],
                      plot_data_df[plot_data_df$mutation_mapped_type == "Ref",],
                      plot_data_df[plot_data_df$mutation_mapped_type == "Var",])
# make colors -------------------------------------------------------------
## make color palette for different read types
colors_mutation_mapped_type <- c("#E31A1C", "#33A02C", "grey70")
names(colors_mutation_mapped_type) <- c("Var", "Ref", "NA")
colors_idents <- c(Polychrome::palette36.colors(n = 36), Polychrome::dark.colors(n = 5))
names(colors_idents) <- as.character(0:40)
# make plots --------------------------------------------------------------

## ggplot
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df[plot_data_df$mutation_mapped_type != "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=mutation_mapped_type), alpha = 0.5, size = 0.3)
p <- p + geom_point_rast(data = plot_data_df[plot_data_df$mutation_mapped_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=mutation_mapped_type), alpha = 0.8, size = 0.8)
p <- p + scale_color_manual(values = colors_mutation_mapped_type)
# p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
p <- p +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p <- p + theme(legend.position = "top")
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())
## save as png
file2write <- paste0(dir_out, "mutationmappedtype_on_umap.", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

#