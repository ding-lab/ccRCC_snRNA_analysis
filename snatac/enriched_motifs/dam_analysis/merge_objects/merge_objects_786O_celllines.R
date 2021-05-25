# Yige Wu @WashU May 2021
## source activate ccrcc_snrna
## reference: https://satijalab.org/signac/articles/merging.html

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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
## library additional libaries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

# set up running parameters -----------------------------------------------
###some parallelization-solution from the tutorial:
future::plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM

# input individual RDS objects into a list ---------------------------------------------------
samples <- c("BAP1_786O", "Control_786O")
paths_rds <- c("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/4.Cell_lines/BAP1-20210412-786-O_processed_atac.chromvar.rds",
               "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/4.Cell_lines/control-20210412-786-O_processed_atac.chromvar.rds")
atac=vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
  atac[[i]]=paths_rds[i]
  DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
  atac[[i]][['RNA']]<-NULL
  atac[[i]][['peaks']]<-NULL
  print (paste(i,samples[i],sep=' '))
  atac[[i]]$Piece_ID=samples[i]
}


# Creating a common peak set ----------------------------------------------
## If the peaks were identified independently in each experiment then they will likely not overlap perfectly. We can merge peaks from all the datasets to create a common peak set, and quantify this peak set in each experiment prior to merging the objects.
#####To obtain the best results - use ALL peaks!
combined.peaks <- Signac::UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- GenomicRanges::width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
####Now using MACS2-peak calling:
## for testing only: subsampling
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)


# Quantify peaks in each dataset ------------------------------------------
## We can now create a matrix of peaks x cell for each sample using the FeatureMatrix function. This function is parallelized using the future package. See the parallelization vignette for more information about using future.
#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.
matrix.counts=vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
  matrix.counts[[i]] <- Signac::FeatureMatrix(
    fragments = Signac::Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
  ) 
}

# create the seurat objects -----------------------------------------------
## We will now use the quantified matrices to create a Seurat object for each dataset, storing the Fragment object for each dataset in the assay.
for (i in 1:length(samples)){
  atac[[i]][['peaksinters']] <- Signac::CreateChromatinAssay(counts = matrix.counts[[i]],
                                                             fragments=Signac::Fragments(atac[[i]]@assays$X500peaksMACS2))
  atac[[i]]$dataset=samples[i]
  DefaultAssay(atac[[i]])<-'peaksinters'
  ###remove other assay
  atac[[i]][['X500peaksMACS2']]<-NULL
}

# Merge objects -----------------------------------------------------------
## Now that the objects each contain an assay with the same set of features, we can use the standard merge function to merge the objects. This will also merge all the fragment objects so that we retain the fragment information for each cell in the final merged object.
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
## save output
file2write <- paste0(dir_out, "786O_CellLines.Merged.", run_id, ".RDS")
saveRDS(object = combined, file = file2write, compress = T)


