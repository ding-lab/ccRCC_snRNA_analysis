# Yige Wu @WashU Dec 2020

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
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_slotdata_bycellgroup_w_epithelialcelltypes_byaliquot_on_katmai/20201218.v1/avgexp.SCT.bycellgroup_w_epithelialcelltypes.byaliquot.20201218.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify parameters ------------------------------------------------------
genes2filter <- "PBRM1"
easyids_snatac <- c("C3N-00733-T1", "C3L-00610-T1", "C3L-01313-T1", "C3L-00416-T2", "C3L-01287-T1", "C3L-00917-T1", "C3L-00088-T1", "C3N-01200-T1", "C3L-00088-T2", "C3L-00079-T1", "C3L-00448-T1",
                    "C3L-00088-N", "C3N-01200-N")
aliquots_snatac <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac]

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycellgroup_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_bycellgroup_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cellgroup = str_split_fixed(string = id_bycellgroup_byaliquot, pattern = "_", n = 2)[,2])
plot_data_long_df$easy_id <- mapvalues(x = plot_data_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_long_df$sample_type <- mapvalues(x = plot_data_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Sample_Type))
plot_data_long_df <- plot_data_long_df %>%
  filter(easy_id %in% easyids_snatac) %>%
  filter((cellgroup == "Tumor.cells" & sample_type == "Tumor") | (cellgroup == "Proximal.tubule" & sample_type == "Normal")) %>%
  mutate(sample_type_text = ifelse(sample_type == "Normal", "NAT", "Tumor"))
  arrange(desc(value))
plot_data_long_df$easy_id_ordered <- factor(x = plot_data_long_df$easy_id, levels = plot_data_long_df$easy_id)

# make barplot ------------------------------------------------------------
p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = easy_id_ordered, y = value))
p <- p + facet_grid(.~sample_type_text, scales = "free_x", shrink = T, space = "free_x")
p <- p + theme_classic()
p <- p + ylab(label = "Average snRNA expression")
p <- p + theme(strip.background = element_rect(fill = NA),
               panel.spacing = unit(0, "lines"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, genes2filter, "_snRNA_expression_across_snATAC_samples.", "pdf")
pdf(file2write, width = 3, height = 3, useDingbats = F)
print(p)
dev.off()

