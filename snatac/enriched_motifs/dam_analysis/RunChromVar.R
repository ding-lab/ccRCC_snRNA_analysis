####Important!!! --to limit number of cores:
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(30)) #number of cores to use - otherwise it crushes

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

#atac=readRDS('21_ccRCC_snATAC.recalledPeaks.v2.20210424.rds')
DefaultAssay(atac)='peaksMACS2'
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(atac),
  pwm = pfm,
  genome = 'BSgenome.Hsapiens.UCSC.hg38',
  use.counts = FALSE
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
atac <- SetAssayData(
  object = atac,
  assay = 'peaksMACS2',
  slot = 'motifs',
  new.data = motif
)

atac[["peaksMACS2"]]

atac <- RegionStats(object = atac, genome = BSgenome.Hsapiens.UCSC.hg38)

atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atac,"26_ccRCC_snATAC.recalledPeaks.chromvar.v3.20210503.rds")