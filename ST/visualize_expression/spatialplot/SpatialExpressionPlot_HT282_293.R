## THE START FROM SCARTCH VERSION
library(tidyverse)
library(Seurat)
library(patchwork)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
library(ggpubr)

## if package is not installed, it will be installed here
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# Load
ffpe_path = '~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/Spatial_transcriptomics/output_FFPE'
stobjs = list()
stobjs$ht282 = Load10X_Spatial(str_glue('{ffpe_path}/HT282N1/H3/HT282N1-S1H3Fs4U1Bp1/outs'))
stobjs$ht293 = Load10X_Spatial(str_glue('{ffpe_path}/HT293N1/H3/HT293N1-S1H3Fs1U1Bp1/outs'))

#######################################
### Quick QC [Optional]
## Note - FFPE no MT genes. MT zero
stobjs = map(stobjs, function(obj){
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)
  })

#######################################
# SCTransform
# https://satijalab.org/seurat/articles/spatial_vignette.html
stobjs = map(stobjs, SCTransform, assay = "Spatial", verbose = T)

#######################################
# plot HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION genes --------------------------------------------------------------------
genes_check <- str_split(string = "CXCL6/IL6/SPP1/CXCL1/NNMT/THBS2/PRRX1/CDH6/MYL9/CXCL8/VCAM1/SERPINE1/INHBA/LUM/CD44/TGFBI/NT5E/ANPEP/AREG/FMOD/TGM2/TFPI2/THBS1/ITGB3/COL6A3/TNFRSF11B/CAPG/IGFBP3/LGALS1/SPARC/GREM1/COL16A1/SPOCK1/ADAM12/IL32/MXRA5/COL5A2/OXTR/MYLK/FGF2/MMP14/DKK1/CCN2/BDNF/PDGFRB/LOXL2/TPM1/LAMC2/FN1/LAMA2/SFRP4/MATN2/GLIPR1/SAT1/ECM1/TNFRSF12A/PLOD2/ITGB1/VEGFC/TPM2/SDC4/CCN1/DAB2/IL15/PTX3/RGS4/VIM/QSOX1/MMP2/MEST/ITGAV/CD59/COL4A1/VCAN/EFEMP2/FSTL3/LRP1/BASP1/TAGLN/COL7A1/IGFBP4/DPYSL3/FBN1/PMEPA1/FAS/PLAUR/CALU/FSTL1/ITGA5/SDC1/COL6A2/ITGB5/COL4A2/EMP3/LOX/DST/CAP2/ENO2/MATN3/CALD1/LAMA1/SERPINH1/BMP1/ACTA2/COPA/PVR/TPM4/LAMA3/GEM/PLOD1/PPIB", pattern = "\\/")[[1]]
dir_out_tmp <- paste0(dir_out, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "/")
dir.create(dir_out_tmp)
for (gene_tmp in genes_check) {
  genes_plot_tmp <- c("CA9", "CP", gene_tmp)
  # file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p1)
  # dev.off()
  
  # file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p2)
  # dev.off()
  
  file2write <- paste0(dir_out_tmp, "HT293_283.", gene_tmp, "_CP", ".png")
  png(file2write, width = 750, height = 700, res = 150)
  p <- ggarrange(p1, p2, nrow = 2)
  print(p)
  dev.off()
}

# plot HALLMARK_INTERFERON_GAMMA_RESPONSE genes --------------------------------------------------------------------
genes_check <- str_split(string = "IL6/SAMD9L/CASP1/PSMB8/XAF1/TNFAIP6/FPR1/VCAM1/PARP14/PARP12/TNFSF10/CFH/OAS2/GBP4/SLAMF7/IFI44/CCL5/IFIH1/APOL6/CASP4/NMI/IL7/RTP4/TNFAIP2/CCL2/PSMB9/C1S/IRF5/BST2/MX1/MVP/IFITM3/UBE2L6/C1R/LGALS3BP/IFIT2/ICAM1/VAMP5/MYD88/MT2A/STAT1/IFI27/SP110/IFIT3/IFI44L/NLRC5/HLA-G/VAMP8/CD274/CD74/DHX58/SOCS3/CD40/PTPN6/SOD2/IL15/IFI35/SSPN/HIF1A/STAT4/IRF1/OASL/PTGS2/DDX58/EPSTI1/DDX60/PSMB10/HLA-A/FAS/CASP8/TXNIP/B2M/PML/TRIM25/CDKN1A/HLA-B/PFKP/TRIM21/LATS2/ISG15/SRI/IRF2/IRF7/NAMPT/ZNFX1/NFKBIA/PLSCR1/IL4R/PSME1/STAT3/RIPK2/TAP1/LAP3/PSME2", pattern = "\\/")[[1]]
dir_out_tmp <- paste0(dir_out, "HALLMARK_INTERFERON_GAMMA_RESPONSE", "/")
dir.create(dir_out_tmp)
for (gene_tmp in genes_check) {
  genes_plot_tmp <- c("CA9", "CP", gene_tmp)
  # file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p1)
  # dev.off()
  
  # file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p2)
  # dev.off()
  
  file2write <- paste0(dir_out_tmp, "HT293_283.", gene_tmp, "_CP", ".png")
  png(file2write, width = 750, height = 700, res = 150)
  p <- ggarrange(p1, p2, nrow = 2)
  print(p)
  dev.off()
}

# plot HALLMARK_TNFA_SIGNALING_VIA_NFKB genes --------------------------------------------------------------------
genes_check <- str_split(string = "CXCL6/IL6/CSF2/CXCL1/TLR2/TNFAIP6/IL18/BIRC3/SERPINE1/INHBA/CD44/CXCL2/PLAU/CCL5/AREG/IL7R/CEBPD/IFIH1/PLK2/PLEK/IL1A/TNFAIP2/G0S2/CCL2/PTPRE/PHLDA1/IL1B/KYNU/IER3/CLCF1/LAMB3/IFIT2/ICAM1/MSC/LIF/ABCA1/SERPINB8/SAT1/DRAM1/SQSTM1/CCL20/SOCS3/SDC4/FOSL1/CCN1/SIK1/SOD2/PTX3/FOSL2/CSF1/SLC2A3/ZC3H12A/TNFSF9/TNIP1/IRF1/SMAD3/CEBPB/PTGS2/HBEGF/DDX58/STAT5A/TGIF1/BCL3/PMEPA1/PFKFB3/RELB/CXCL3/PLAUR/SGK1/IL6ST/TNFRSF9/CFLAR/CDKN1A/MAFF/TIPARP/NAMPT/TUBB2A/TNFAIP8/TNIP2/NFKBIA/LITAF/EHD1/BIRC2/NFE2L2/TANK/KLF6/NFKBIE/RIPK2/TAP1/NFIL3/GEM/MCL1/B4GALT1/NFAT5", pattern = "\\/")[[1]]
dir_out_tmp <- paste0(dir_out, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "/")
dir.create(dir_out_tmp)
for (gene_tmp in genes_check) {
  genes_plot_tmp <- c("CA9", "CP", gene_tmp)
  # file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p1)
  # dev.off()
  
  # file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p2)
  # dev.off()
  
  file2write <- paste0(dir_out_tmp, "HT293_283.", gene_tmp, "_CP", ".png")
  png(file2write, width = 750, height = 700, res = 150)
  p <- ggarrange(p1, p2, nrow = 2)
  print(p)
  dev.off()
}

# plot HALLMARK_INFLAMMATORY_RESPONSE genes --------------------------------------------------------------------
genes_check <- str_split(string = "EREG/CXCL6/IL6/TLR2/AQP9/TNFAIP6/FPR1/IL18/OSMR/CXCL8/SERPINE1/TNFSF10/INHBA/TNFRSF1B/AXL/CX3CL1/CCL5/MSR1/IL7R/ITGB3/IL1A/NMI/RTP4/ITGB8/TLR3/CCL2/PTPRE/IL1B/NLRP3/BST2/DCBLD2/ADGRE1/MMP14/SELL/ICAM1/LIF/ABCA1/CD82/AHR/IL1R1/CD70/KCNJ2/RNF144B/TLR1/CCL20/HRH1/GPR132/MET/CD40/IL15/CSF1/IRAK2/HIF1A/TNFSF9/IRF1/HBEGF/C3AR1/BEST1/PLAUR/ITGA5/TNFRSF9/EBI3/CDKN1A/EMP3/SRI/IRF7/NAMPT/RHOG/NFKBIA/LPAR1/IL4R/PVR/KLF6/RIPK2/SGMS2/ADRM1", pattern = "\\/")[[1]]
dir_out_tmp <- paste0(dir_out, "HALLMARK_INFLAMMATORY_RESPONSE", "/")
dir.create(dir_out_tmp)
for (gene_tmp in genes_check) {
  genes_plot_tmp <- c("CA9", "CP", gene_tmp)
  # file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p1)
  # dev.off()
  
  # file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p2)
  # dev.off()
  
  file2write <- paste0(dir_out_tmp, "HT293_283.", gene_tmp, "_CP", ".png")
  png(file2write, width = 750, height = 700, res = 150)
  p <- ggarrange(p1, p2, nrow = 2)
  print(p)
  dev.off()
}

# plot NABA_MATRISOME_ASSOCIATED genes --------------------------------------------------------------------
genes_check <- str_split(string = "EREG/CXCL6/IL6/CSF2/CXCL5/CLEC2B/HRNR/CXCL1/PI3/SERPINB7/COLEC10/CLEC4E/FLG/PF4V1/FLG2/IL18/S100A16/S100A6/CXCL8/SERPINE1/ADAMTS16/TNFSF10/INHBA/ANXA13/MMP7/BTC/SERPINA1/CXCL2/SERPINB9/S100A3/TNFSF14/PLAU/CX3CL1/CCL5/CTSZ/CTSS/GDF15/AREG/TGM2/ADAM28/FGF1/BMP5/IL1A/LGALS1/IL7/GREM1/TGFA/CCL2/CCL28/ADAMTS12/IL1B/ANGPTL4/P3H2/ADAM12/CLCF1/FGF2/MMP14/ANXA1/S100A2/BDNF/LOXL2/TGFB2/SLPI/LIF/SFRP4/HGF/MMP24/SERPINB8/ADAMTSL1/LGALS3/P4HA3/ANXA2/MMP19/GDF5/PAPPA2/EPGN/PLOD2/CCL20/VEGFC/FGFBP1/PLXNB3/S100A13/IL34/P4HTM/SDC4/TIMP2/SEMA3C/IL15/S100A11/CTSD/CTSB/CSF1/NRG1/CD109/MUC12/LGALS8/MMP2/ANXA4/TNFSF9/FSTL3/S100A4/P4HA1/HBEGF/PAPPA/SEMA4B/CXCL3/P4HA2/ANXA5/FSTL1/SDC1/ADAMTS3/EBI3/ELFN2/EGLN1/ANXA3/LOX/CTSO/FGF5/PDGFC/CTSC/ADAM9/ANXA6/ANXA7/SERPINH1/BMP1/ADAM10/ANXA11/S100A10/ADAM15/SERPINB6/CTSA/PLOD1/ADAM17", pattern = "\\/")[[1]]
dir_out_tmp <- paste0(dir_out, "NABA_MATRISOME_ASSOCIATED", "/")
dir.create(dir_out_tmp)
for (gene_tmp in genes_check) {
  genes_plot_tmp <- c("CA9", "CP", gene_tmp)
  # file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p1)
  # dev.off()
  
  # file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
  # png(file2write, width = 750, height = 350, res = 150)
  p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
  # print(p2)
  # dev.off()
  
  file2write <- paste0(dir_out_tmp, "HT293_283.", gene_tmp, "_CP", ".png")
  png(file2write, width = 750, height = 700, res = 150)
  p <- ggarrange(p1, p2, nrow = 2)
  print(p)
  dev.off()
}


# test --------------------------------------------------------------------
genes_plot_tmp <- c("CA9", "CP", "COL4A1", "OSMR", "OSM", "TGM2", "VEGFA")
# file2write <- paste0(dir_out_tmp, "HT283.", gene_tmp, "_CP", ".png")
# png(file2write, width = 750, height = 350, res = 150)
p1 <- SpatialPlot(stobjs$ht282, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
# print(p1)
# dev.off()

# file2write <- paste0(dir_out_tmp, "HT293.", gene_tmp, "_CP", ".png")
# png(file2write, width = 750, height = 350, res = 150)
p2 <- SpatialPlot(stobjs$ht293, features = genes_plot_tmp, stroke= NA, image.alpha=0, ncol = length(genes_plot_tmp))
# print(p2)
# dev.off()

file2write <- paste0(dir_out, "HT293_283.", paste0(genes_plot_tmp, collapse = "_"), ".png")
png(file2write, width = 1600, height = 600, res = 150)
p <- ggarrange(p1, p2, nrow = 2)
print(p)
dev.off()

# Spatial Plot
genes_check = c("LAMA1", "LAMA3", "COL5A1", "COL8A1", "COL28A1", "CD70", "CD27", "OSM", "OSMR", "IL6ST")
genes_check = c("COL5A2", "COL6A3", "COL16A1", "COL8A1", "COL5A2", "COL28A1", "COL4A1", "COL7A1", "COL6A2", "COL4A2", "COL18A1", "COL17A1", "COL5A1",
                "COL12A1", "COL1A1", "COL4A3", "COL27A1", "COL6A1", "COL10A1", "COL4A4", "COL4A4",
                "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4")

file2write <- paste0(dir_out, "HT283.collagens.", "png")
png(file2write, width = 2000, height = 3000, res = 150)
SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

file2write <- paste0(dir_out, "HT293.collagens.", "png")
png(file2write, width = 2000, height = 3000, res = 150)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

genes_check <- c("LAMB3", "LAMC2", "LAMA2", "LAMB1", "LAMA1", "LAMA3", "LAMC1", "LAMA4", "LAMA5", "LAMB2")
file2write <- paste0(dir_out, "HT283.", paste0(genes_check, collapse = "_"), ".png")
png(file2write, width = 1400, height = 800, res = 150)
SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

file2write <- paste0(dir_out, "HT293.", paste0(genes_check, collapse = "_"), ".png")
png(file2write, width = 1400, height = 800, res = 150)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

genes_check <- str_split(string = "EREG/CXCL6/IL6/CSF2/CXCL5/CLEC2B/HRNR/CXCL1/PI3/SERPINB7/COLEC10/CLEC4E/FLG/PF4V1/FLG2/IL18/S100A16/S100A6/CXCL8/SERPINE1/ADAMTS16/TNFSF10/INHBA/ANXA13/MMP7/BTC/SERPINA1/CXCL2/SERPINB9/S100A3/TNFSF14/PLAU/CX3CL1/CCL5/CTSZ/CTSS/GDF15/AREG/TGM2/ADAM28/FGF1/BMP5/IL1A/LGALS1/IL7/GREM1/TGFA/CCL2/CCL28/ADAMTS12/IL1B/ANGPTL4/P3H2/ADAM12/CLCF1/FGF2/MMP14/ANXA1/S100A2/BDNF/LOXL2/TGFB2/SLPI/LIF/SFRP4/HGF/MMP24/SERPINB8/ADAMTSL1/LGALS3/P4HA3/ANXA2/MMP19/GDF5/PAPPA2/EPGN/PLOD2/CCL20/VEGFC/FGFBP1/PLXNB3/S100A13/IL34/P4HTM/SDC4/TIMP2/SEMA3C/IL15/S100A11/CTSD/CTSB/CSF1/NRG1/CD109/MUC12/LGALS8/MMP2/ANXA4/TNFSF9/FSTL3/S100A4/P4HA1/HBEGF/PAPPA/SEMA4B/CXCL3/P4HA2/ANXA5/FSTL1/SDC1/ADAMTS3/EBI3/ELFN2/EGLN1/ANXA3/LOX/CTSO/FGF5/PDGFC/CTSC/ADAM9/ANXA6/ANXA7/SERPINH1/BMP1/ADAM10/ANXA11/S100A10/ADAM15/SERPINB6/CTSA/PLOD1/ADAM17", pattern = "\\/")[[1]]
file2write <- paste0(dir_out, "HT283.matrisome.", "png")
png(file2write, width = 1500, height = 5000, res = 150)
SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

file2write <- paste0(dir_out, "HT293.matrisome.", "png")
png(file2write, width = 1500, height = 5000, res = 150)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

genes_check = c("OSM", "OSMR", "IL6ST")
genes_check = c("CP", "VCAM1", "VIM", "ICAM1")
genes_check = c("CP", "ANGPTL4", "VEGFC", "VEGFA")
genes_check = c("CA9", "CP", "COL4A2", "COL4A1", "LOX")
genes_check = c("CA9", "CP", "COL4A1")
genes_check = c("CA9", "CP", "COL4A1", "LAMA5")

file2write <- paste0(dir_out, "HT283.", paste0(genes_check, collapse = "_"), ".png")
png(file2write, width = 250*length(genes_check), height = 350, res = 150)
SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0, ncol = length(genes_check))
dev.off()

file2write <- paste0(dir_out, "HT293.", paste0(genes_check, collapse = "_"), ".png")
png(file2write, width = 250*length(genes_check), height = 350, res = 150)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = length(genes_check))
dev.off()






