# Yige Wu @WashU Sep 2019
## make bubbleplot for non-immune clusters like Fig 2A: https://www.nature.com/articles/s41467-019-10861-2

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# functions ---------------------------------------------------------------
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    result <- FetchData(object[["RNA"]], genes)
    result <- colMeans(result  > 0)*100
    result <- data.frame(gene_symbol = names(result), perc_exp = result)
    return(result)
  } else if (group.by == "seurat_clusters") {
    clusters <- unique(object@meta.data$seurat_clusters)
    
  } else {        
    stop("not specified group!")
  }
}


# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input the integrated non-immune clusters --------------------------------
renal.int.nonimmune.obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/recluster/process_nonimmune_clusters/20191003.v1/NonImmune_Integrated.20191003.v1.RDS")

# input the DEG table -----------------------------------------------------
renal.markers <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/recluster/process_nonimmune_clusters/20191003.v1/NonImmune.DEGs.Pos.CellTypeMarkerOnly.20191003.v1.txt", data.table = F)

# # input the gene2cell type annotation file --------------------------------
# gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191003.v1.tsv", data.table = F)

# get the list of cell type marker genes differentially expressed --------------------
renal.markers$gene %>% unique() %>% length()
## 77 genes
genes2plot <- renal.markers$gene %>% unique()

# calculate percent of cells expressing marker genes for each cluster --------------
DefaultAssay(renal.int.nonimmune.obj) <- "RNA"
p <- DotPlot(object = renal.int.nonimmune.obj, features = genes2plot, col.min = 0)
p <- p  + RotatedAxis()
file2write <- paste0(dir_out, biomarker_type_tmp, "_Marker_Expression.", run_id, ".pdf")
pdf(file2write, 
    width = 12,
    height = 16)
p



