# Yige Wu @WashU Jun 2020

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
## input barcode2manualclusterid
barcode2manualcluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv")
which(is.na(barcode2manualcluster_df$Id_TumorManualCluster))
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)

# make data frame to plot------------------------
## filter out barcode without subcluster assignment
plot_data_df <- barcode2manualcluster_df %>%
  filter(!is.na(Id_TumorManualCluster)) %>%
  select(orig.ident, Id_TumorManualCluster) %>%
  unique() %>%
  select(orig.ident) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  rename(Id_Aliquot = '.')
## add reable id
plot_data_df$Id_Aliquot_WU <- mapvalues(x = plot_data_df$Id_Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## order by frequency
plot_data_df <- plot_data_df %>%
  arrange(Freq)
plot_data_df$plot_x <- factor(x = plot_data_df$Id_Aliquot_WU, levels = rev(plot_data_df$Id_Aliquot_WU))

# plot wide version --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, mapping = aes(x = plot_x, y = Freq), stat = "identity")
p <- p + ylab("Number of Tumor Subclusters")
p <- p + ggtitle("Number of Tumor Subclusters Per Sample")
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank())
p <- p + theme(axis.text.y = element_text(size = 15),
               axis.text.x = element_text(angle = 90, vjust = 0.1))
p
pdf2write <- paste0(dir_out, "NumberOfManualTumorSubclusters.wide", ".pdf")
pdf(file = pdf2write, width = 8, height = 5)
print(p)
dev.off()


# plot long version --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, mapping = aes(x = plot_x, y = Freq), stat = "identity")
p <- p + coord_flip()
p <- p + ylab("Number of Tumor Subclusters")
p <- p + scale_y_continuous(position = "right")
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 12))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
p <- p + theme(axis.line = element_line(colour = "black"), axis.ticks.y = element_blank())
p
pdf2write <- paste0(dir_out, "NumberOfManualTumorSubclusters.long",".pdf")
pdf(file = pdf2write, width = 4, height = 6, useDingbats = F)
print(p)
dev.off()

