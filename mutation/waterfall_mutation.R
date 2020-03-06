# Yige Wu @WashU March 2020
## show the bulk mutation and 10Xmapping result in a walterfall plot

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## library
library(ggplot2)
# input dependencies ------------------------------------------------------
## input bulk mutation data
bulk_mutation_df <- loadMaf()

## input 10xmapping result
snRNA_mutation_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)

## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)

# make the data frame for plotting ----------------------------------------
## get the cases to show
cases2show <- unique(meta_tab$Case[meta_tab$Aliquot.snRNA %in% srat_paths$Aliquot])
cases2show
## reformat bulk mutation data
### create a column for case id to filter later
bulk_mutation_df <- bulk_mutation_df %>%
  mutate(Case = stringr::str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_T", n = 2)[,1])
### filter by case
bulk_mutation_df <- bulk_mutation_df %>%
  filter(Case %in% cases2show)
### filter by SMGs
bulk_mutation_df <- bulk_mutation_df %>%
  filter(Hugo_Symbol %in% ccRCC_SMGs)
### select specific columns
bulk_mutation_df <- bulk_mutation_df %>%
  select(Case, Hugo_Symbol, HGVSp_Short)

## reformat the snRNA mutation mapping
### get the aliquot is for the original tumor piece
snRNA_aliquot_ids <- meta_tab$Aliquot.snRNA[meta_tab$Is_discovery_set == T & meta_tab$Case %in% cases2show & meta_tab$Sample_Type == "Tumor"]
snRNA_aliquot_ids
### filter by the original tumor piece
snRNA_mutation_df <- snRNA_mutation_df %>%
  filter(aliquot %in% snRNA_aliquot_ids)
### collapse the number of cells with detected variant read to each aliquot
snRNA_mutation_cell_count <- snRNA_mutation_df %>%
  filter(allele_type == "Var") %>%
  group_by(aliquot, gene_symbol, mutation) %>%
  summarise(num_cells = length(unique(barcode)))
### filter by SMGs
snRNA_mutation_cell_count_filtered <- snRNA_mutation_cell_count %>%
  filter(gene_symbol %in% ccRCC_SMGs)
### map aliquot to case id
snRNA_mutation_cell_count_filtered$Case <- plyr::mapvalues(x = snRNA_mutation_cell_count_filtered$aliquot, from = meta_tab$Aliquot.snRNA, to = meta_tab$Case)

## merge bulk mutation with snRNA mutation mapping
plot_data_df <- merge(bulk_mutation_df, snRNA_mutation_cell_count_filtered, by.x = c("Case", "Hugo_Symbol"), by.y = c("Case", "gene_symbol"), all.x = T)
## order x and y axis
plot_data_df <- plot_data_df %>%
  mutate(Var1 = Case) %>%
  mutate(Var2 = Hugo_Symbol) %>%
  mutate(Value = 1)
### get the matrix from the long data frame
plot_data_mat <- reshape2::dcast(data = plot_data_df, formula = Var2 ~ Var1 , fill = 0, value.var = "Value")
rownames(plot_data_mat) <- plot_data_mat$Var2
plot_data_mat <- plot_data_mat %>%
  select(-Var2)
plot_data_mat[plot_data_mat!=0] <- 1
#define a graph that represented as adjacency matrix with matrix A
library(igraph)
g <- graph.incidence(plot_data_mat, weighted = TRUE)
# cluster wit Louvain algorithm
lou <- cluster_louvain(g)
df.lou <- data.frame(lou$names,lou$membership)
#the same longData than earlier
plot_data_df <- left_join(plot_data_df, df.lou, by=c("Var1"="lou.names"))
colnames(plot_data_df)[ncol(plot_data_df)] <- "Var1_clust"
plot_data_df$Var2 <- as.factor(plot_data_df$Var2)
plot_data_df <- left_join(plot_data_df, df.lou, by=c("Var2"="lou.names"))
colnames(plot_data_df)[ncol(plot_data_df)] <- "Var2_clust"
plot_data_df$Var1 <- factor(plot_data_df$Var1, levels=unique(arrange(plot_data_df, Var1_clust)[,1]))
plot_data_df$Var2 <- factor(plot_data_df$Var2, levels=unique(arrange(plot_data_df, Var2_clust)[,2]))
plot_data_df$Var2 <- factor(plot_data_df$Var2, levels=rev(ccRCC_SMGs))

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_tile(data = plot_data_df, mapping = aes(x = Var1, y = Var2), fill = "lightblue")
p <- p + geom_point(data = plot_data_df[!is.na(plot_data_df$num_cells),], mapping = aes(x = Var1, y = Var2))
p <- p + xlab("Case ID")
p <- p + ylab("Mutated Genes")
p <- p + theme_grid()
p <- p + theme(axis.text.x = element_text(size = 10, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p

# save plot ---------------------------------------------------------------
png(filename = paste0(dir_out, "Grid_Bulk_Mutation_and_10XMapping", ".", run_id, ".", ".png"), width = 1000, height = 300, res = 150)
print(p)
dev.off()

