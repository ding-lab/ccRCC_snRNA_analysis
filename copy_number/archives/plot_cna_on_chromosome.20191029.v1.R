# Yige Wu @WashU Oct 2019
## plotting copy number region from subcluster mode of infercnv


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
packages = c(
  "karyoploteR"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# create output directory -------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input infercnv subcluster cnv region results ----------------------------
snRNA_aliquot_id_tmp <- "CPT0001260013"
cnv_region_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.20191021.v1/CPT0001260013/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat", data.table = F)
cnv_region_tab <- cnv_region_tab %>%
  mutate(cluster_id = str_split_fixed(string = cell_group_name, pattern = "\\.", n = 6)[,1])


# create output directory by aliquot --------------------------------------
dir_out_aliquot <- paste0(dir_out, snRNA_aliquot_id_tmp, "/")
dir.create(dir_out_aliquot)

# plot by subcluster on each chromosome--------------------------------------------------------------------
for (chr_tmp in c("chr11")) {
  dir_out_chr <- paste0(dir_out_aliquot, chr_tmp, "/")
  dir.create(dir_out_chr)
  
  for (subcluster_id_tmp in unique(cnv_region_tab$cell_group_name)) {
    cnv_region2plot <- cnv_region_tab %>%
      filter(cell_group_name == subcluster_id_tmp) %>%
      filter(chr == chr_tmp)
    file2write <- paste(dir_out_chr, snRNA_aliquot_id_tmp, "_Subcluster", subcluster_id_tmp, "_", chr_tmp, "_CNA.", run_id, ".png", sep="")
    png(file2write, width = 600, height = 500, res = 150)
    kp <- plotKaryotype(chromosomes=c(chr_tmp), main = paste0("Subcluster ", subcluster_id_tmp, "\nsnRNA-Inferred Copy Number Ratios At ", chr_tmp))
    kpAddBaseNumbers(kp)
    kpDataBackground(kp, data.panel = 1, r0=-0.05, r1=1.05)
    kpAxis(kp, data.panel = 1, ymin = 0, ymax = 2, r0 = 0, r1 = 2/2, numticks = 5)
    for (i in 1:nrow(cnv_region2plot)) {
      start_tmp <- cnv_region2plot$start[i]
      end_tmp <- cnv_region2plot$end[i]
      state_tmp <- cnv_region2plot$state[i]
      kpLines(kp, data.panel = 1, chr = chr_tmp, x = c(start_tmp, end_tmp), y = c(state_tmp, state_tmp)/2, pch=".", cex=2)
    }
    kpAbline(kp, v=c(23016671), col = "red")
    kpText(kp, chr=chr_tmp, x=23016671, y=1/2, col="red", labels="t(3,11)", cex=0.7)
    dev.off()
  }
}

for (chr_tmp in c("chr5")) {
  dir_out_chr <- paste0(dir_out_aliquot, chr_tmp, "/")
  dir.create(dir_out_chr)
  
  for (subcluster_id_tmp in unique(cnv_region_tab$cell_group_name)) {
    cnv_region2plot <- cnv_region_tab %>%
      filter(cell_group_name == subcluster_id_tmp) %>%
      filter(chr == chr_tmp)
    file2write <- paste(dir_out_chr, snRNA_aliquot_id_tmp, "_Subcluster", subcluster_id_tmp, "_", chr_tmp, "_CNA.", run_id, ".png", sep="")
    png(file2write, width = 600, height = 500, res = 150)
    kp <- plotKaryotype(chromosomes=c(chr_tmp), main = paste0("Subcluster ", subcluster_id_tmp, "\nsnRNA-Inferred Copy Number Ratios At ", chr_tmp))
    kpAddBaseNumbers(kp)
    kpDataBackground(kp, data.panel = 1, r0=-0.05, r1=1.05)
    kpAxis(kp, data.panel = 1, ymin = 0, ymax = 2, r0 = 0, r1 = 2/2, numticks = 5)
    for (i in 1:nrow(cnv_region2plot)) {
      start_tmp <- cnv_region2plot$start[i]
      end_tmp <- cnv_region2plot$end[i]
      state_tmp <- cnv_region2plot$state[i]
      kpLines(kp, data.panel = 1, chr = chr_tmp, x = c(start_tmp, end_tmp), y = c(state_tmp, state_tmp)/2, pch=".", cex=2)
    }
    kpAbline(kp, v=c(87092375, 80708698, 86181276), col = "red")
    kpText(kp, chr=chr_tmp, x=80708698, y=1/2, col="red", labels="t(3,5)", cex=0.7)
    # kpText(kp, chr=chr_tmp, x=86181276, y=1/2, col="red", labels="t(3,5)-B", cex=0.7)
    # kpText(kp, chr=chr_tmp, x=87092375, y=1/2, col="red", labels="t(3,5)-C", cex=0.7)
    dev.off()
  }
}



for (chr_tmp in c("chr3")) {
  dir_out_chr <- paste0(dir_out_aliquot, chr_tmp, "/")
  dir.create(dir_out_chr)
  
  for (subcluster_id_tmp in unique(cnv_region_tab$cell_group_name)) {
    cnv_region2plot <- cnv_region_tab %>%
      filter(cell_group_name == subcluster_id_tmp) %>%
      filter(chr == chr_tmp)
    file2write <- paste(dir_out_chr, snRNA_aliquot_id_tmp, "_Subcluster", subcluster_id_tmp, "_", chr_tmp, "_CNA.", run_id, ".png", sep="")
    png(file2write, width = 600, height = 500, res = 150)
    kp <- plotKaryotype(chromosomes=c(chr_tmp), main = paste0("Subcluster ", subcluster_id_tmp, "\nsnRNA-Inferred Copy Number Ratios At ", chr_tmp))
    kpAddBaseNumbers(kp)
    kpDataBackground(kp, data.panel = 1, r0=-0.05, r1=1.05)
    kpAxis(kp, data.panel = 1, ymin = 0, ymax = 2, r0 = 0, r1 = 2/2, numticks = 5)
    for (i in 1:nrow(cnv_region2plot)) {
      start_tmp <- cnv_region2plot$start[i]
      end_tmp <- cnv_region2plot$end[i]
      state_tmp <- cnv_region2plot$state[i]
      kpLines(kp, data.panel = 1, chr = chr_tmp, x = c(start_tmp, end_tmp), y = c(state_tmp, state_tmp)/2, pch=".", cex=2)
    }
    kpAbline(kp, v=c(81005946, 114492369, 114578911, 114846179), col = "red")
    kpText(kp, chr="chr3", x=81005946, y=1/2, col="red", labels="t(3,11)", cex=0.7)
    kpText(kp, chr="chr3", x=114492369, y=1/2, col="red", labels="t(3,5)", cex=0.7)
    dev.off()
  }
  
}
