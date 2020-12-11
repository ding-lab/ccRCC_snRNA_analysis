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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20201207.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20201207.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# make data frame for plotting --------------------------------------------
## add aliquot.wu
cnv_3state_count_aliquots$aliquot.wu <- mapvalues(x = cnv_3state_count_aliquots$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_df <- cnv_3state_count_aliquots
plot_data_df$case <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

## add cytoband and expected cna type
plot_data_df$gene_cytoband <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
plot_data_df$gene_expected_state <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
plot_data_df <- plot_data_df %>%
  mutate(id_aliquot_cluster = tumor_subcluster) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,2]) %>%
  filter(name_tumorsubcluster != "CNA")
  
## get the data with only expected CNV state
plot_data_cc_df <- plot_data_df %>%
  filter(gene_expected_state == cna_3state) 

## filter genes
genes_filtered <- unique(plot_data_cc_df$gene_symbol[plot_data_cc_df$Fraction > 0.1])
plot_data_cc_df <- plot_data_cc_df %>%
  filter(gene_symbol %in% genes_filtered) %>%
  mutate(Fraction_Range = ifelse(Fraction < 0.5, "<=50%", ">50%")) %>%
  mutate(aliquot.wu = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,1]) %>%
  mutate(Data_detected = T) %>%
  select(id_aliquot_cluster, gene_symbol, Fraction, name_tumorsubcluster, cna_3state, Fraction_Range, Data_detected, aliquot.wu, gene_cytoband)
  
## make data frame for the NA data
all_data_df <- plot_data_df %>%
  filter(gene_symbol %in% genes_filtered) %>%
  select(id_aliquot_cluster, gene_symbol, gene_cytoband, aliquot.wu) %>%
  unique()
plot_data_na_df <- all_data_df %>%
  select(id_aliquot_cluster, gene_symbol) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 0) %>%
  mutate(Data_detected = F) %>%
  mutate(Fraction = 1) %>%
  mutate(aliquot.wu = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,1]) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,2]) %>%
  mutate(cna_3state = "NA") %>%
  mutate(Fraction_Range = "NA") %>%
  select(id_aliquot_cluster, gene_symbol, Fraction, name_tumorsubcluster, cna_3state, Fraction_Range, Data_detected, aliquot.wu)
plot_data_na_df$gene_cytoband <- mapvalues(x = plot_data_na_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
## merge
plot_data_df <- rbind(plot_data_cc_df, plot_data_na_df)

# order x axis facet ------------------------------------------------------
## count the number of subclusters per aliquot
count_subclusters_df <- cnv_3state_count_aliquots %>%
  select(aliquot.wu, tumor_subcluster) %>%
  unique() %>%
  select(aliquot.wu) %>%
  table() %>%
  as.data.frame() %>%
  rename(aliquot.wu = '.') %>%
  rename(count_subclusters = Freq) %>%
  arrange(count_subclusters)
levels_aliquot.wu <- count_subclusters_df$aliquot.wu
## order aliquot id by transformating into factor
plot_data_df$aliquot.wu <- factor(plot_data_df$aliquot.wu, levels = levels_aliquot.wu)

# make data frame for highlighted gene and sample -------------------------
hl_data_df <- plot_data_df %>%
  group_by(aliquot.wu, gene_symbol) %>%
  summarize(Fraction_maxdiff = (max(Fraction) - min(Fraction)), x1 = id_aliquot_cluster[which.min(Fraction)], x2 = id_aliquot_cluster[which.max(Fraction)]) %>%
  filter(Fraction_maxdiff >= 0.75)
plot_data_comb_df <- merge(plot_data_df, hl_data_df, by.x = c("id_aliquot_cluster", "gene_symbol", "aliquot.wu"), by.y = c("x1", "gene_symbol", "aliquot.wu"), all.x = T)
table(plot_data_comb_df$gene_cytoband)
plot_data_comb_df$gene_cytoband <- factor(x = plot_data_comb_df$gene_cytoband, levels = unique(knowncnvgenes_df$Cytoband))

# plot expected CNVs wide--------------------------------------------------------------------
p <- ggplot(data = plot_data_comb_df)
p <- p + geom_point(mapping = aes(y = gene_symbol, 
                                  x = id_aliquot_cluster, 
                                  size = Fraction, 
                                  fill = cna_3state, 
                                  color = Fraction_Range,
                                  shape = Data_detected),
                    alpha = 0.8)
p <- p + scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 4))
p <- p + scale_fill_manual(values = c("Gain" = "red", "Loss" = "blue", "NA" = "#000000"))
p <- p + scale_color_manual(values = c(">50%" = "#000000", "<=50%" = "#FFFFFF", "NA" = "#000000"))
p <- p + geom_segment(aes(x = id_aliquot_cluster, y = gene_symbol, xend = x2, yend = gene_symbol), arrow = arrow(length = unit(0.2,"cm")))
p <- p + facet_grid(gene_cytoband~aliquot.wu,
                    scales = "free", space = "free", shrink = T)
p <- p + scale_x_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5),
               axis.text.y = element_text(size = 12, face = "bold"))
p <- p + theme(panel.spacing = unit(0, "lines"))
p <- p + theme(strip.text.y = element_text(angle = 0),
               strip.text.x = element_text(angle = 90))

pdf(paste0(dir_out, "Fraction_of_Tumorcells_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".pdf"), 
    width = 23, height = 6)
print(p)
dev.off()

png(file = paste0(dir_out, "Fraction_of_Tumorcells_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".png"), 
    width = 3500, height = 1000, res = 150)
print(p)
dev.off()

# plot expected CNVs long--------------------------------------------------------------------
p <- ggplot(data = plot_data_comb_df)
p <- p + geom_point(mapping = aes(x = gene_symbol, 
                                  y = id_aliquot_cluster, 
                                  size = Fraction, 
                                  fill = cna_3state, 
                                  color = Fraction_Range,
                                  shape = Data_detected),
                    alpha = 0.8)
p <- p + scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 4))
p <- p + scale_fill_manual(values = c("Gain" = "red", "Loss" = "blue", "NA" = "#000000"))
p <- p + scale_color_manual(values = c(">50%" = "#000000", "<=50%" = "#FFFFFF", "NA" = "#000000"))
# p <- p + geom_segment(aes(y = id_aliquot_cluster, x = gene_symbol, xend = x2, yend = gene_symbol), arrow = arrow(length = unit(0.2,"cm")))
p <- p + facet_grid(aliquot.wu~gene_cytoband,
                    scales = "free", space = "free", shrink = T) # , switch = "y"
p <- p + scale_y_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + theme_bw()
p <- p + theme(axis.title = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0, hjust = 0.5),
               axis.text.y = element_text(size = 12, face = "bold"))
p <- p + theme(panel.spacing = unit(0, "lines"))
p <- p + theme(strip.text.x = element_text(angle = 90),
               strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
               strip.background = element_rect(color="black", fill="white"))
p <- p + theme(legend.position = "left")

pdf(paste0(dir_out, "Fraction_of_Tumorcells_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".", "long", ".pdf"), 
    width = 8, height = 20)
print(p)
dev.off()

png(file = paste0(dir_out, "Fraction_of_Tumorcells_with_Expected_CNA_by_Gene", ".", "By_Tumor_Subcluster", ".", "long",  ".png"), 
    width = 1000, height = 3000, res = 150)
print(p)
dev.off()

