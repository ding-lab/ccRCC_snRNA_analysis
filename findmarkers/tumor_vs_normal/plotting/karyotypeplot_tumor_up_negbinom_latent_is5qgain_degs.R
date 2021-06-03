# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(biomaRt)
library(karyoploteR)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_negbinom_is5qgain_all_ccRCC_vs_pt_on_katmai/20210602.v1/negbinom.latent.Is.5q.Gain.logfc.threshold0.25.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")

# make plot data ----------------------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(avg_log2FC > 0) %>%
  filter(p_val_adj < 0.05)
genes2convert <- unique(deg_filtered_df$genesymbol_deg)
genesymbol2region_df <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', "start_position", "end_position", "strand", "band"), 
                              filters = 'hgnc_symbol', 
                              values = genes2convert, 
                              mart = ensembl)
genesymbol2region_df <- genesymbol2region_df %>%
  arrange(hgnc_symbol, desc(band))
plotdata_df <- genesymbol2region_df[!duplicated(genesymbol2region_df$hgnc_symbol),]
plotdata_df <- plotdata_df %>%
  mutate(gene_region = paste0("chr", chromosome_name,"-", start_position, "-", end_position))
regions=Signac::StringToGRanges(plotdata_df$gene_region, sep = c("-", "-"))

#  plot----------------------------------------------------------
## plot
file2write <- paste0(dir_out, "Up.png")
png(file2write, width = 1000, height = 800, res = 150)
kp <- plotKaryotype(genome="hg38")
kpPlotRegions(kp, regions, col="red")
dev.off()