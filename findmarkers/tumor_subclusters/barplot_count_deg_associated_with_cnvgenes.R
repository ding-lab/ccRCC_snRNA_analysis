# Yige Wu @WashU May 2020

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

# input dependencies -------------------------------------------------------
## load CNV fraction in tumor cells
cnv_frac_df <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200512.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200512.v1.tsv", data.table = F)
## input count of degs
count_deg_by_ppi_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_tumormanualsubcluster_deg_by_ppi/20200601.v1/tumormanualsubcluster.deg_count.cnv_genes_interactome.20200601.v1.tsv")

# get all the cluster names --------------------------------------------
## add aliquot.wu
cnv_frac_df$aliquot.wu <- mapvalues(x = cnv_frac_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add cytoband and expected cna type
cnv_frac_df$gene_cytoband <- mapvalues(x = cnv_frac_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
cnv_frac_df$gene_expected_state <- mapvalues(x = cnv_frac_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
cnv_frac_df <- cnv_frac_df %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1)))

# create plot data --------------------------------------------------------
plot_data_df <- count_deg_by_ppi_df %>%
  mutate(text = ifelse(name_cluster == "C3L-00088-T2_C2", paste0("(", cna_gene_symbol, ")"), NA)) %>%
  mutate(Freq.plot = ifelse(deg_direction == "down", -Freq, Freq)) %>%
  mutate(aliquot.wu = str_split_fixed(string = name_cluster, pattern = "_", n = 2)[,1]) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = name_cluster, pattern = "_", n = 2)[,2])

test <- count_deg_by_ppi_df %>%
  group_by(name_cluster, deg_direction) %>%
  top_n(1)
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = name_cluster, y = Freq.plot, color = deg_direction, size = sqrt(Freq)), shape = 16, alpha = 0.7)
p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = name_cluster, y = Freq.plot, label = text), angle = 90, size = 3)
p <- p + scale_color_manual(values = c("up" = "red", "down" = "blue"))
p <- p + facet_grid(.~aliquot.wu,
                    scales = "free", space = "free", shrink = T)
p <- p + scale_x_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 12, face = "bold"))
p <- p + theme(panel.spacing = unit(0, "lines"))
p <- p + theme(strip.text.y = element_text(angle = 0),
               strip.text.x = element_text(angle = 90))
p
