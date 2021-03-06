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
## library additional libaries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# set up running parameters -----------------------------------------------
###some parallelization-solution from the tutorial:
future::plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM
print("Finished setting running parameters")

# input individual RDS objects into a list ---------------------------------------------------
samples <- c("BAP1_786O", "Control_786O")
paths_rds <- c("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/4.Cell_lines/BAP1-20210412-786-O_processed_atac.chromvar.rds",
               "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/4.Cell_lines/control-20210412-786-O_processed_atac.chromvar.rds")
atac=vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
  atac[[i]]=readRDS(file = paths_rds[i])
  Seurat::DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
  atac[[i]][['RNA']]<-NULL
  atac[[i]][['peaks']]<-NULL
  print (paste(i,samples[i],sep=' '))
  atac[[i]]$Piece_ID=samples[i]
}
print("Finished creating the atac list")

# Creating a common peak set ----------------------------------------------
## If the peaks were identified independently in each experiment then they will likely not overlap perfectly. We can merge peaks from all the datasets to create a common peak set, and quantify this peak set in each experiment prior to merging the objects.
#####To obtain the best results - use ALL peaks!
combined.peaks <- Signac::UnifyPeaks(object.list = atac, mode = "reduce")
print("Finished UnifyPeaks")
peakwidths <- GenomicRanges::width(combined.peaks)
print("Finished width")
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
print("Finished filtering combined.peaks")
# combined.peaks
####Now using MACS2-peak calling:
## 5K peaks are not the ones in the final object this step is just to make initial merge (to have fragment files together) but you will use the peaks from the first script
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)
# peaks.use <- combined.peaks

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
print("Finished Quantify peaks")

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
print("Finished creating seurat objects")

# Merge objects -----------------------------------------------------------
## Now that the objects each contain an assay with the same set of features, we can use the standard merge function to merge the objects. This will also merge all the fragment objects so that we retain the fragment information for each cell in the final merged object.
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
print("Finished merge")

# Create a gene activity matrix -------------------------------------------
gene.activities <- GeneActivity(combined)
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

## save output
file2write <- paste0(dir_out, "786O_CellLines.Merged.", run_id, ".RDS")
saveRDS(object = combined, file = file2write, compress = T)
print("Finished saveRDS")


