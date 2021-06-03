# Yige Wu @WashU Jun 2021
## source activate ccrcc_snrna

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
## library additional libaries
library(Signac)
library(Seurat)
library(future)
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
plan("multiprocess", workers = 10)
options(future.globals.maxSize= 891289600)

# input dependencies ------------------------------------------------------
atac=readRDS('Resources/Analysis_Results/snatac/merge_objects/add_gene_activity_to_786O_celllines_merged_katmai/20210528.v1/786O_CellLines.Merged.20210528.v1.RDS')
recentered_final <- read.table(file = "./Resources/Analysis_Results/snatac/overlap_peaks/overlap_peaks_786O_celllines/20210601.v1/recentered_final.filtered.20210601.v1.tsv", sep='\t',header=TRUE)

# create featurematrix with the new set of peaks --------------------------
recentered_p=Signac::StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- Signac::FeatureMatrix(
  fragments = Signac::Fragments(atac@assays$peaksinters),
  features = recentered_p,
  sep = c("-","-"),
  cells = colnames(atac)
)
atac[['peaksMACS2']] <- Signac::CreateChromatinAssay(counts = matrix.counts,
                                             fragments=Signac::Fragments(atac@assays$peaksinters))
SeuratObject::DefaultAssay(atac)<-'peaksMACS2'

atac[['X500peaksMACS2']]<-NULL
atac[['peaksinters']]<-NULL
print(paste0("Finished creating peaksMACS2 assay!"))

# clustering --------------------------------------------------------------
atac <- RunTFIDF(atac)
print(paste0("Finished RunTFIDF!"))

atac <- FindTopFeatures(atac, min.cutoff = 'q0')
print(paste0("Finished FindTopFeatures!"))

atac <- RunSVD(atac)
print(paste0("Finished RunSVD!"))

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:50)
print(paste0("Finished RunUMAP!"))

atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:50)
print(paste0("Finished FindNeighbors!"))

atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
print(paste0("Finished FindNeighbors!"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "786O_CellLines.PeakOverlapped.20210528.v1.RDS")
saveRDS(object = atac, file = file2write, compress = T)


