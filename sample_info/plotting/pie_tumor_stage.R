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
clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_case_clinical_data/20200717.v1/snRNA_ccRCC_Clinicl_Table.20200717.v1.tsv")

# make plot data frame ----------------------------------------------------
plot_data_df <- clinical_df %>%
  select(Tumor_Stage_Pathological) %>%
  table() %>%
  as.data.frame() %>%
  rename(Tumor_Stage_Pathological = '.') %>%
  mutate(Fraction = (Freq/nrow(clinical_df))*100) %>%
  arrange(desc(Fraction)) %>%
  mutate(ypos = cumsum(Fraction)- 0.5*Fraction) %>%
  mutate(text_label = paste0("n = ", Freq))
# mutate(text_label = paste0(round(x = Fraction, digits = 1), "%"))
plot_data_df$Tumor_Stage_Pathological <- factor(x = plot_data_df$Tumor_Stage_Pathological, levels = rev(plot_data_df$Tumor_Stage_Pathological))


# make colors -------------------------------------------------------------
colors_stage <- RColorBrewer::brewer.pal(n = 8, name = "YlOrRd")[c(2, 4, 6, 8)]
swatch(colors_stage)
names(colors_stage) <- c("Stage I", "Stage II", "Stage III", "Stage IV")
# plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_df, aes(x="", y = Fraction, fill=Tumor_Stage_Pathological))
p <- p + geom_bar(stat="identity", width=1, color = "white")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values = colors_stage)
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
