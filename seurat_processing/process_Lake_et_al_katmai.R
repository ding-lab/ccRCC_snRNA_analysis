
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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input dependencies ------------------------------------------------------
input <- fread(data.table = F, input  = "../RCC_Literature_Review/Single_Cell/2019_sn-Seq_HumanKidney_NatComm_Jain/GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv.gz")
dim(input)
cluster_anno_df <- fread(data.table = F, input = "../RCC_Literature_Review/Single_Cell/2019_sn-Seq_HumanKidney_NatComm_Jain/GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv.gz")

# process to seurat object ------------------------------------------------
input_mat <- as.matrix(input[,-1])
rownames(input_mat) <- input$V1
rm(input)
srat = CreateSeuratObject(counts = input_mat, project= "Lake et al.",min.cells = 0)
rm(input_mat)
##
srat@meta.data$cluster_Lake <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 3)[,1]
table(srat@meta.data$cluster_Lake)
##
cluster_anno_df <- cluster_anno_df %>%
  mutate(Cluster_Name = paste0("C", Cluster))
srat@meta.data$Abbn <- mapvalues(x = srat@meta.data$cluster_Lake, from = cluster_anno_df$Cluster_Name, to = as.vector(cluster_anno_df$Abbn))
table(srat@meta.data$Abbn)
##
srat@meta.data$library <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 3)[,2]
##
mito.genes <- grep(pattern = "^MT-", x = rownames(srat), value=F)
percent.mito <- Matrix::colSums((GetAssayData(srat,slot="counts"))[mito.genes, ])/Matrix::colSums(GetAssayData(srat,slot="counts"))
srat$percent.mito<-percent.mito
##
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA","percent.mito", "library"), return.only.var.genes = F)
srat <- RunPCA(srat, npcs = 30, verbose = FALSE)
cat("Finished RUNPCA!\n")
srat <- RunUMAP(srat, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat <- FindClusters(srat, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
file2write <- paste0(dir_out, "Lake_etal.Seurat.RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("Finished saving the output!\n")

