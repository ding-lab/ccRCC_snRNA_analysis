#####author: Nadezhda V. Terekhanova
#####2020-09-26: v.3.0 - this version of the script was adaapted for the latest Signac v.1.0.0
#####2020-07-01##################
#####v.2.0 - "auto", now just need to provide the sample-ids. Script will create the list of Seurat objects, and will store the provided sample ids in the $dataset (meta-data).
#####Output - Seurat object of the merged datasets and list of samples' ids used.

###############
#IMPORTANT: do prior to running the script `export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp` - otherwise it crushes (probably because not enough space to write everything)
###############

library(Signac)
library(Seurat)
library(GenomeInfoDb)
###library(EnsDb.Hsapiens.v75) ###Now using newer version of the annotation:
library(ggplot2)
library(RColorBrewer)

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
library(GenomicRanges)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM

samples=c('CPT0000880001', 'CPT0000890002', 'CPT0001260013', 'CPT0010100001', 'CPT0010160004',
 'CPT0023690004', 'CPT0025110004', 'CPT0025880013', 'CPT0075130004', 'CPT0075170013', 'CPT0079410004',
'CPT0086820004', 'CPT0001180011', 'CPT0012280004', 'CPT0014450005','CPT0000870003','CPT0001220012',
'CPT0012550012','CPT0063630004','CPT0075720013','CPT0078510004','CPT0001500003', 'CPT0001540013',
 'CPT0019130004', 'CPT0065690004', 'CPT0086350004')


tab=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'snRNA_processed_Yige/meta_data.20200505.v1.tsv',sep=''),sep='\t',header=TRUE)
tab=tab[tab$Aliquot.snRNA %in% samples,]
samples=tab$Aliquot.snRNA
piece_ids=tab$Aliquot.snRNA.WU

atac=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    atac[[i]]=readRDS(paste("../../1.Create_rds/out/",samples[i],"/",samples[i],"_processed_atac.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['RNA']]<-NULL
    atac[[i]][['peaks']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}

#####To obtain the best results - use ALL peaks!

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#peaks.use=combined.peaks
####Now using MACS2-peak calling:
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#For testing purposes only:
#peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
    ) 
}


for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],
fragments=Fragments(atac[[i]]@assays$X500peaksMACS2))
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
###remove other assay
atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
saveRDS(combined, '26_ccRCC_snATAC.20210503.rds')





