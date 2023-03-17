# Yige Wu @WashU Jan 2022
## plot mutation mapping type on integration UMAP

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
table(plot_data_df$mutation_mapped_type)

### order the data frame so that cells mapped with variant allele will come on top
plot_data_var_df <- plot_data_df[plot_data_df$mutation_mapped_type == "Var",]
plot_data_var_df$genesymbol_mapped <- mapvalues(x = plot_data_var_df$id_barcode, from = barcodes_mut_df$id_barcode, to = as.vector(barcodes_mut_df$gene_symbol))
plot_data_ref_df <- plot_data_df[plot_data_df$mutation_mapped_type == "Ref",]
plot_data_ref_df$genesymbol_mapped <- mapvalues(x = plot_data_ref_df$id_barcode, from = barcodes_ref_df$id_barcode, to = as.vector(barcodes_ref_df$gene_symbol))
plot_data_nocov_df <- plot_data_df[plot_data_df$mutation_mapped_type == "NA",]
plot_data_nocov_df$genesymbol_mapped <- NA
plot_data_df <- rbind(plot_data_nocov_df,
                      plot_data_ref_df,
                      plot_data_var_df)
plot_data_df <- plot_data_df %>%
  mutate(mutation_mapped_type2 = ifelse(mutation_mapped_type == "Var",
                                        ifelse(genesymbol_mapped %in% ccRCC_SMGs, "Var_SMG", "Var_nonSMG"),
                                        ifelse(mutation_mapped_type == "Ref", "Ref", "NA")))
table(plot_data_df$mutation_mapped_type2)
plot_data_df <- plot_data_df %>%
  arrange(factor(x = mutation_mapped_type2, levels = rev(c("Var_SMG", "Var_nonSMG", "Ref", "NA"))))
plot_data_df$ident <- as.character(plot_data_df$ident)

# make colors -------------------------------------------------------------
## make color palette for different read types
colors_mutation_mapped_type2 <- c("#E31A1C", "#FF7F00", "#33A02C", "grey70")
names(colors_mutation_mapped_type2) <- c("Var_SMG", "Var_nonSMG", "Ref", "NA")

colors_mutation_mapped_type <- c("#E31A1C", "#33A02C", "grey70")
names(colors_mutation_mapped_type) <- c("Var", "Ref", "NA")
colors_idents <- c(Polychrome::palette36.colors(n = 36), Polychrome::dark.colors(n = 5))
names(colors_idents) <- as.character(0:40)

# make plots --------------------------------------------------------------
## ggplot
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df[plot_data_df$mutation_mapped_type2 != "Var_SMG",], mapping = aes(x = UMAP_1, y = UMAP_2, color=mutation_mapped_type2), alpha = 0.5, size = 0.3)
p <- p + geom_point_rast(data = plot_data_df[plot_data_df$mutation_mapped_type2 == "Var_SMG",], mapping = aes(x = UMAP_1, y = UMAP_2, color=mutation_mapped_type2), alpha = 1, size = 1.5)
p <- p + scale_color_manual(values = colors_mutation_mapped_type2)
p <- p +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 1, label.theme = element_text(size = 18)))
p <- p + theme(legend.position="bottom", aspect.ratio=1)

## save as pdf
file2write <- paste0(dir_out, "mutationmappedtype_on_umap.", "pdf")
pdf(file = file2write, width = 6, height = 6.5, useDingbats = F)
print(p)
dev.off()

## save source data
source_data_df <- plot_data_df %>%
  select(UMAP_1, UMAP_2, mutation_mapped_type2)
write.table(x = source_data_df, file = "~/Desktop/SF1c.SourceData.tsv", quote = F, sep = "\t", row.names = F)
