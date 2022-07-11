# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/ora_msigdb_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/ora_msigdb_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
enrich_df <- fread(data.table = F, input = "~/Downloads/ORA_Results.tsv")

# make plot data ----------------------------------------------------------
count_geneset_df <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(Description, easy_id) %>%
  unique() %>%
  dplyr::select(Description) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Description = ".") %>%
  arrange(desc(Freq))
plotdata_df <- count_geneset_df %>%
  mutate(GeneSet_Name = gsub(x = Description, pattern = "HALLMARK_", replacement = "")) %>%
  head(38)
plotdata_df$GeneSet_Name[plotdata_df$GeneSet_Name == "EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT"
plotdata_df$GeneSet_Name[plotdata_df$GeneSet_Name == "REACTIVE_OXYGEN_SPECIES_PATHWAY"] <- "REACTIVE_OXYGEN_SPECIES"
plotdata_df$GeneSet_Name <- factor(x = plotdata_df$GeneSet_Name, levels = plotdata_df$GeneSet_Name)

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = GeneSet_Name, y = Freq), stat = "identity", fill = RColorBrewer::brewer.pal(name = "Accent", n = 4)[1])
p <- p + geom_text(data = plotdata_df, mapping = aes(x = GeneSet_Name, y = 0.5, label = GeneSet_Name), angle = 90, hjust = "bottom", size = 3)
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())
p <- p + theme(axis.text.y = element_text(size = 15), axis.line.y = element_blank())
p <- p + ylab("No. tumors with over-represented \n DEGs by gene set")
p
file2write <- paste0(dir_out, "Count.ORA.P.adjust.0.05.png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "Count.ORA.P.adjust.0.05.pdf")
pdf(file2write, width = 8, height = 2.75, useDingbats = F)
print(p)
dev.off()
