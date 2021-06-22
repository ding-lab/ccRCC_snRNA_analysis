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
                          Is_differentially_expressed = c("TRUE", "FALSE", "TRUE", "FALSE"))


# plot  -------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = peak2gene_type, y = No_peaks, fill = Is_differentially_expressed), stat = "identity")
p <- p + scale_fill_manual(values = c("TRUE" = "yellow", "FALSE" = "grey40"))
p <- p + theme_classic()
p <- p + theme(axis.title.x = element_blank())
# p <- p + theme(legend.position = "")
p
