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
## input gene-to-pathway
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pathway/unite_PBRM1_BAP1_vs_NonMutants_DAP_ORA/20210628.v1/PBRM1_BAP1_vs_NonMutants.DAPGene2TopPathway.20210628.v1.tsv")
## input pbrm1 daps
daps_pbrm1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/annotate_pbrm1_vs_nonmutant_daps/20210625.v1/PBRM1_DAP2Gene.EnhancerPromoter.20210625.v1.tsv")
daps_bap1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_peaks/annotate_BAP1_vs_NonMutant_daps/20210625.v1/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20210625.v1.tsv")
## input 

# make plot data ----------------------------------------------------------
colnames_merge <- intersect(colnames(daps_pbrm1_df), colnames(daps_bap1_df))
plotdata_source_df <- rbind(daps_pbrm1_df[,colnames_merge] %>%
                       mutate(DAP_group = paste0("PBRM1_", DAP_direction)),
                     daps_bap1_df[,colnames_merge] %>%
                       mutate(DAP_group = paste0("BAP1_", DAP_direction)))
plotdata_source_df <- merge(x = plotdata_source_df, y = gene2pathway_df, by.x = c("Gene"), by.y = c("GeneSymbol"), all.x = T)  

plotdata_df <- plotdata_source_df %>%
  filter(!is.na(GeneSet_Name)) %>%
  group_by(DAP_group) %>%
  select(GeneSet_Name) %>%
  table() %>%
  as.data.frame()

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq))
p <- p + scale_size_area()
p
file2write <- paste0(dir_out, "bubble.plain.png")
png(file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()
