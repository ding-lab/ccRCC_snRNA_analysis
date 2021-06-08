# Yige Wu @WashU Jun 2021

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
## input annotated table
deg2region_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_pt_snRNA_degs_to_chr_regions/20210602.v1/Tumor_vs_PT.snRNA_DEGs.Chromosome_Regions.20210602.v1.tsv")
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210608.v1/Consistent.Tumor_vs_PT_DEGs.CNVcorrected.20210608.v1.tsv")

#  plot----------------------------------------------------------
## make plot data
plotdata_df <- deg2region_df %>%
  filter(Tumor_vs_PT == "Up") %>%
  filter(hgnc_symbol %in% degs_df$genesymbol_deg) %>%
  select(chromosome_name) %>%
  table() %>%
  as.data.frame() %>%
  rename(chromosome_name = '.')
plotdata_df$chromosome_name <- factor(x = plotdata_df$chromosome_name, levels = c(as.character(1:22), "X", "MT"))
## plot
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = chromosome_name, y = Freq), stat = "identity")
p <- p + theme_classic(base_size = 12)
p <- p + ylim(c(0, 125))
p <- p + ggtitle(paste0("Significant DEGS across >= 15 Tumorcell-vs-PT comparisons"), subtitle = 'CNV corrected')
file2write <- paste0(dir_out, "Up_DEGs.bychr.CNVcorrected.png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()

deg2region_df %>%
  filter(Tumor_vs_PT == "Up") %>%
  filter(!(hgnc_symbol %in% degs_df$genesymbol_deg)) %>%
  select(chromosome_name) %>%
  table()
