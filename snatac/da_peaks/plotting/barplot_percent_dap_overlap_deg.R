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

# make plot data ----------------------------------------------------------
plotdata_df <- data.frame(peak2gene_type = c("Promoter peak", "Promoter peak", "Enhancer peak", "Enhancer peak"),
                          No_peaks = c(73, 269, 144, 212),
                          Peak_for_differentially_expressed_gene = c("TRUE", "FALSE", "TRUE", "FALSE"))

# plot  -------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = peak2gene_type, y = No_peaks, fill = Peak_for_differentially_expressed_gene), stat = "identity")
# p <- p + geom_text(data = plotdata_df, mapping = aes(x = peak2gene_type, y = 0.5, label = peak2gene_type), angle = 90, hjust = "bottom", size = 3)
p <- p + scale_fill_manual(values = c("TRUE" = brewer.pal(n = 12, name = "Set3")[5], "FALSE" = brewer.pal(n = 12, name = "Set3")[12]))
p <- p + theme_classic(base_size = 14)
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1),
               legend.title = element_text(size = 12))
# p <- p + theme(legend.position = "")
p
file2write <- paste0(dir_out, "No_peaks_de.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()