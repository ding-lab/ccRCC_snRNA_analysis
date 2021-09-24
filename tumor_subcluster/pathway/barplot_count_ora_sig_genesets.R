# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
count_geneset_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/count_ora_sig_genesets_across_samples/20210413.v1/Count.ORA.P.adjust.0.05.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- count_geneset_df %>%
  filter(!(GeneSet_Name %in% c("Retinoblastoma Gene in Cancer", "Gastric Cancer Network 1", "Gastric Cancer Network 2", 
                               "Integrated Cancer Pathway", "Integrated Breast Cancer Pathway", "Pathways Affected in Adenoid Cystic Carcinoma"))) %>%
  head(20)
plotdata_df$GeneSet_Name[plotdata_df$GeneSet_Name == "DNA IR-Double Strand Breaks (DSBs) and cellular response via ATM"] <- "DNA IR-Double Strand Breaks (DSBs) &\n cellular response via ATM"
plotdata_df$GeneSet_Name <- factor(x = plotdata_df$GeneSet_Name, levels = plotdata_df$GeneSet_Name)
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = GeneSet_Name, y = Freq), stat = "identity", fill = "lightgreen")
p <- p + geom_text(data = plotdata_df, mapping = aes(x = GeneSet_Name, y = 0.5, label = GeneSet_Name), angle = 90, hjust = "bottom", size = 3)
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())
p <- p + theme(axis.text.y = element_text(size = 15), axis.line.y = element_blank())
p <- p + ylab("No. tumors with over-represented gene set\nin DEGs of their tumor subclusters")
p
file2write <- paste0(dir_out, "Count.ORA.P.adjust.0.05.png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()

