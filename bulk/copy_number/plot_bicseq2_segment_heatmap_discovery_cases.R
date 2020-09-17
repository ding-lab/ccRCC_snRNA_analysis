# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(copynumber)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependecies -------------------------------------------------------



# get file paths to process -----------------------------------------------
dir_tmp <- paste0("./Resources/Bulk_Processed_Data/WGS_CNV_Somatic/BICSEQ2/")
cnv_file_names <- list.files(path = dir_tmp, recursive = T)
cnv_file_names
cnv_file_names <- cnv_file_names[grepl(pattern = "cnv", x = cnv_file_names)]
cnv_file_names
cnv_tab4heatmap <- NULL
for (cnv_file_name in cnv_file_names) {
  cnv_file_path <- paste0(dir_tmp, cnv_file_name)
  if (file.info(cnv_file_path)["size"] > 0) {
    cnv_file <- fread(input = cnv_file_path, data.table = F)
    case <- strsplit(x = cnv_file_name, split = "\\.")[[1]][1]
    chr_col <- str_split_fixed(string = cnv_file$chrom, pattern = "chr", 2)[,2]
    chr_col[chr_col %in% c("X", "Y")] <- "23"
    chr_col <- as.integer(as.vector(chr_col))
    cnv_tab_tmp <- data.frame(sampleID = case, chrom = chr_col, arm = "NA", 
                              start.pos = cnv_file$start, end.pos = cnv_file$end, 
                              n.probes = cnv_file$binNum, mean = cnv_file$log2.copyRatio)
    
    cnv_tab4heatmap <- rbind(cnv_tab4heatmap, cnv_tab_tmp)
  }
}
cnv_tab4heatmap$sampleID <- factor(cnv_tab4heatmap$sampleID, levels = sort(as.vector(unique(cnv_tab4heatmap$sampleID))))
cnv_tab4heatmap <- cnv_tab4heatmap[order(cnv_tab4heatmap$sampleID),]
## plot heatmap
filename <- paste0(dir_out, "heatmap.png")
png(filename, width = 2000, height = 2000, res = 150)
# plotHeatmap(segments=cnv_tab4heatmap,upper.lim=c(0.2))
plotHeatmap(segments=cnv_tab4heatmap,upper.lim=c(0.2), colors=c("dodgerblue","white","red"))
dev.off()

filename <- paste0(dir_out, "heatmap.pdf")
pdf(filename, width = 10, height = 13)
plotHeatmap(segments=cnv_tab4heatmap,upper.lim=c(0.2), colors=c("dodgerblue","white","red"))
dev.off()

