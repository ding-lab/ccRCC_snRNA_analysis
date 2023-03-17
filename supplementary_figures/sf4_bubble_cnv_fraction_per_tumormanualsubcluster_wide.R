# Yige Wu @WashU May 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
# idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv", data.table = F)
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster_rm_doublets/20210806.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20210806.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input clusters to plotf
clusters_selected_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/plotting/heatmap/heatmap_hallmark38_scores_by_manual_tumorcluster_sct_data_scaled_selected/20220706.v1/intrapatient_tumorclusters.selected.20220706.v1.tsv")

# preprocess --------------------------------------------------------------
clusters_selected_df <- clusters_selected_df %>%
  mutate(id_aliquot_cluster = gsub(pattern = "\\.", replacement = "-", x = cluster_name))

# make data frame for plotting --------------------------------------------
## add aliquot.wu
cnv_3state_count_aliquots$aliquot.wu <- mapvalues(x = cnv_3state_count_aliquots$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_df <- cnv_3state_count_aliquots
plot_data_df$case <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
plot_data_df <- plot_data_df %>%
  filter(case != "C3L-00359")
## add cytoband and expected cna type
plot_data_df$gene_cytoband <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
plot_data_df$gene_expected_state <- mapvalues(x = plot_data_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
plot_data_df <- plot_data_df %>%
  mutate(id_aliquot_cluster = tumor_subcluster) %>%
  filter(id_aliquot_cluster %in% clusters_selected_df$id_aliquot_cluster) %>%
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
count_subclusters_df <- plot_data_cc_df %>%
  select(aliquot.wu, id_aliquot_cluster) %>%
  unique() %>%
  select(aliquot.wu) %>%
  table() %>%
  as.data.frame() %>%
  rename(aliquot.wu = '.') %>%
  rename(count_subclusters = Freq) %>%
  arrange(count_subclusters)
count_subclusters_df$cumsum_count <- cumsum(count_subclusters_df$count_subclusters)
levels_aliquot.wu <- count_subclusters_df$aliquot.wu
## order aliquot id by transformating into factor
plot_data_df$aliquot.wu <- factor(plot_data_df$aliquot.wu, levels = levels_aliquot.wu)
## filter genes again
df_tmp <- dcast(data = plot_data_df, gene_symbol ~ id_aliquot_cluster, value.var = "cna_3state")
df_tmp2 <- df_tmp[,-1]; rownames(df_tmp2) <- df_tmp$gene_symbol
na_count_bygene <- rowSums(!is.na(df_tmp2) & df_tmp2 == "NA")
genes_keep <- names(na_count_bygene)[na_count_bygene <= nrow(clusters_selected_df)/2]
plot_data_comb_df <- plot_data_df %>%
  filter(gene_symbol %in% genes_keep)
table(plot_data_comb_df$gene_cytoband)
## order cytobands
plot_data_comb_df$gene_cytoband <- factor(x = plot_data_comb_df$gene_cytoband, levels = unique(knowncnvgenes_df$Cytoband))
unique(plot_data_comb_df$id_aliquot_cluster)

# make plotting parameters ------------------------------------------------
color_cn_gain <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_cn_loss <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]

# plot expected CNVs wide - 2 rows--------------------------------------------------------------------
fontsize_plot <- 13

mytheme <- theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = fontsize_plot, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = fontsize_plot, face = "italic", color = "black"))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text.y = element_text(size = 12, angle = 0, color = "black"),
        strip.text.x = element_text(size = fontsize_plot, angle = 90, color = "black"),
        strip.background = element_rect(color="black", fill="white"))+
  theme(axis.title.y = element_blank())
library(gridExtra)
samples_p1 <- count_subclusters_df$aliquot.wu[count_subclusters_df$cumsum_count <= nrow(clusters_selected_df)/2]; samples_p1 <- as.vector(samples_p1)
samples_p2 <- count_subclusters_df$aliquot.wu[count_subclusters_df$cumsum_count > nrow(clusters_selected_df)/2]; samples_p1 <- as.vector(samples_p1)
plotdata_df1 <- subset(plot_data_comb_df, aliquot.wu %in% samples_p1)
plotdata_df2 <- subset(plot_data_comb_df, aliquot.wu %in% samples_p2)

p1 <- ggplot(data = plotdata_df1) +
  geom_point(mapping = aes(y = gene_symbol, x = id_aliquot_cluster, size = Fraction, color = cna_3state, shape = Data_detected), alpha = 0.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 15))+
  scale_color_manual(values = c("Gain" = color_cn_gain, "Loss" = color_cn_loss, "NA" = "grey80"))+
  facet_grid(gene_cytoband~aliquot.wu, scales = "free", space = "free", shrink = T )+
  scale_x_discrete(breaks=plotdata_df1$id_aliquot_cluster, labels=plotdata_df1$name_tumorsubcluster)+
  guides(colour = guide_legend(override.aes = list(size=4), title = "CNV state", title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1),
         shape = guide_legend(override.aes = list(size=4),  title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1),
         size = guide_legend(title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1))+
  mytheme

p2 <- ggplot(data = plotdata_df2) +
  geom_point(mapping = aes(y = gene_symbol, x = id_aliquot_cluster, size = Fraction, color = cna_3state, shape = Data_detected), alpha = 0.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 15))+
  scale_color_manual(values = c("Gain" = color_cn_gain, "Loss" = color_cn_loss, "NA" = "grey80"))+
  facet_grid(gene_cytoband~aliquot.wu, scales = "free", space = "free", shrink = T )+
  scale_x_discrete(breaks=plotdata_df2$id_aliquot_cluster, labels=plotdata_df2$name_tumorsubcluster)+
  guides(colour = guide_legend(override.aes = list(size=4), title = "CNV state", title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1),
         shape = guide_legend(override.aes = list(size=4),  title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1),
         size = guide_legend(title.theme = element_text(size = fontsize_plot), label.theme = element_text(size = fontsize_plot), ncol = 1))+
  mytheme


# write output ------------------------------------------------------------
pdf(paste0(dir_out, "2row", ".pdf"), 
    width = 12, height = 13, useDingbats = F)
grid.arrange(p1, p2, nrow = 2)
dev.off()
## write source data
write.table(x = plot_data_comb_df, file = "~/Desktop/SF4.SourceData.tsv", quote = F, sep = "\t", row.names = F)
