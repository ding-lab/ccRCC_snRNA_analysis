# Yige Wu @WashU Sep 2020

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
## input the case to mutation category
mut_cat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_sample_by_pbrm1_bap1_mutation/20200914.v1/PBRM1_BAP1_Mutation_Status_By_Case.20200914.v1.tsv")

# make plot data frame ----------------------------------------------------
plot_data_df <- mut_cat_df %>%
  select(mutation_category) %>%
  table() %>%
  as.data.frame() %>%
  rename(mutation_category = '.') %>%
  mutate(Fraction = (Freq/nrow(mut_cat_df))*100) %>%
  arrange(desc(Fraction)) %>%
  mutate(ypos = cumsum(Fraction)- 0.5*Fraction) %>%
  mutate(text_label = paste0("n = ", Freq))
  # mutate(text_label = paste0(round(x = Fraction, digits = 1), "%"))
plot_data_df$mutation_category <- factor(x = plot_data_df$mutation_category, levels = rev(plot_data_df$mutation_category))

# plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_df, aes(x="", y = Fraction, fill=mutation_category))
p <- p + geom_bar(stat="identity", width=1, color = "white")
p <- p + coord_polar("y", start=0)
# p <- p + scale_fill_manual(values = normal_epithelial_colors)
# p <- p + scale_fill_brewer(palette="Set1")
p <- p + theme_void()
p <- p + geom_text(aes(x= 1.3, y = ypos, label = text_label), color = "black", size=6)
# p <- p + ggtitle(label = paste0("Cell Type Distribution of ", sum(plot_data_df$value), " Normal Epithelial Cells"))
p <- p + theme(legend.position = "top")
p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
p

# write output -------------------------------------------------------------
file2write <- paste0(dir_out, "pie.pdf")
pdf(file2write, height = 6, width = 6)
print(p)
dev.off()
file2write <- paste0(dir_out, "pie.png")
png(file2write, height = 800, width = 900, res = 150)
print(p)
dev.off()
