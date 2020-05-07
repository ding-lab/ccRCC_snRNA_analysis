# Yige Wu @WashU May 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
# devtools::install_github("teunbrand/ggh4x")
library(ggh4x)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200505.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200505.v1.tsv", data.table = F)
table(cnv_3state_count_aliquots$tumor_subcluster)
## input the CNV type per tumor
cnvtype_aliquot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/summarize_cnv_ith_per_tumor/20200505.v1/CNV_Type_Assignment_Per_Tumor.20200505.v1.tsv")
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# make data frame for plotting --------------------------------------------
## add aliquot.wu
cnv_3state_count_aliquots$aliquot.wu <- mapvalues(x = cnv_3state_count_aliquots$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_df <- cnv_3state_count_aliquots
plot_data_df$case <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
plot_data_df$aliquot_cnv_type <- mapvalues(x = plot_data_df$aliquot.wu, from = cnvtype_aliquot_df$aliquot.wu, to = as.vector(cnvtype_aliquot_df$CNV_ITH_Type))
## add cytoband and expected cna type
plot_data_df$gene_cytoband <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
plot_data_df$gene_expected_state <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
## make data frame for the background
all_data_df <- plot_data_df %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1))) %>%
  filter(gene_symbol %in% y_filtered) %>%
  select(id_aliquot_cluster, gene_symbol, gene_cytoband, aliquot.wu) %>%
  unique()
plotbackground_data_df <- all_data_df %>%
  select(id_aliquot_cluster, gene_symbol) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 0) %>%
  mutate(Data_Notdetected = T)

## filter out copy neutral ones
plot_data_df <- plot_data_df %>%
  filter(gene_expected_state == cna_3state) %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1))) %>%
  mutate(Fraction_Range = ifelse(Fraction < 0.5, "<=50%", ">50%")) %>%
  mutate(segment_suffix = str_split_fixed(string = aliquot.wu, pattern = "-", n = 3)[,3]) %>%
  mutate(name_tumorsubcluster = paste0("C", (tumor_subcluster + 1))) %>%
  filter(cna_3state != "Neutral")
  
## filter rows
y_filtered <- unique(plot_data_df$gene_symbol[plot_data_df$Fraction > 0.5])
plot_data_df <- plot_data_df %>%
  filter(gene_symbol %in% y_filtered)

# make colors --------------------------------------------------------------
cnvtype_aliquot_df <- cnvtype_aliquot_df %>%
  mutate(color_cnv_type = ifelse(CNV_ITH_Type == "Intra-segment_Heterogeneous", "yellow", "grey")) 

# order x axis facet ------------------------------------------------------
cnvtype_aliquot_df <- cnvtype_aliquot_df %>%
  arrange(CNV_ITH_Type) 
levels_aliquot.wu <- cnvtype_aliquot_df$aliquot.wu
plot_data_df$aliquot.wu <- factor(plot_data_df$aliquot.wu, levels = levels_aliquot.wu)
# plot expected CNVs--------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(y = gene_symbol, x = id_aliquot_cluster, size = Fraction, fill = cna_3state, color = Fraction_Range),
                    shape = 21, alpha = 0.8)
# p <- p + scale_fill_manual(values = copy_number_colors)
p <- p + scale_color_manual(values = c(">50%" = "#000000", "<=50%" = "#FFFFFF"))
p <- p + facet_grid(gene_cytoband~aliquot.wu,
                    scales = "free", space = "free", shrink = T)
p <- p + scale_x_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 12, face = "bold"))
p <- p + theme(panel.spacing = unit(0, "lines"))
p <- p + theme(strip.text.y = element_text(angle = 0),
               strip.text.x = element_text(angle = 90))

# add color to the facet --------------------------------------------------
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))

for (i in strip_both) {
  ## get text
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  facet_text <- g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$label
  k <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[k]]$gp$fill <- mapvalues(x = facet_text, from = cnvtype_aliquot_df$aliquot.wu, to = as.vector(cnvtype_aliquot_df$color_cnv_type))
}

# p
png(file = paste0(dir_out, "Fraction_of_Tumorcells_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".", run_id, ".png"), 
    width = 2500, height = 1000, res = 150)
# print(p)
grid.draw(g)
dev.off()
