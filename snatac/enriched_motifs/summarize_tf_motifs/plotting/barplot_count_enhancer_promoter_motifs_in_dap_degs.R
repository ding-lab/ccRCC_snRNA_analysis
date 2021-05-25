# Yige Wu @WashU May 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
count_motifs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/count_motif_promoter_enhancer_occurance_in_up_dap_deg/20210511.v2/Count_Motif_in_DAP_DEGs.20210511.v2.tsv")

# make plot data ----------------------------------------------------------
plotorder_df <- count_motifs_df %>%
  group_by(TF_name) %>%
  summarise(Freq.all = sum(Freq)) %>%
  arrange(desc(Freq.all))
plotdata_df <- count_motifs_df %>%
  filter(TF_name %in% head(x = plotorder_df$TF_name, n = 30)) %>%
  mutate(motif_location = DAP_Type.strict)
plotdata_df$x_plot <- factor(x = plotdata_df$TF_name, levels = plotorder_df$TF_name) 
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = x_plot, y = Freq, fill = motif_location), stat = "identity")
# p <- p + geom_text(data = plotdata_df, mapping = aes(x = TF_name, y = 0.5, label = GeneSet_Name), angle = 90, hjust = "bottom", size = 3)
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.2), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())
p <- p + theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 12), axis.line.y = element_blank(), legend.position=c(.9,.85))
p <- p + ylab("No. genes with the motif in\ndifferentially accessible regions")
p
file2write <- paste0(dir_out, "Count.png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "Count.pdf")
pdf(file2write, width = 5, height = 3, useDingbats = F)
print(p)
dev.off()
