# Yige Wu @WashU Sep 2019
## make bubbleplot for non-immune clusters like Fig 2A: https://www.nature.com/articles/s41467-019-10861-2

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input the integrated non-immune clusters --------------------------------
renal.int.nonimmune.obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/recluster/process_nonimmune_clusters/20191003.v1/NonImmune_Integrated.20191003.v1.RDS")

# input the DEG table -----------------------------------------------------
renal.markers <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/recluster/process_nonimmune_clusters/20191003.v1/NonImmune.DEGs.Pos.CellTypeMarkerOnly.20191003.v1.txt", data.table = F)

# input the gene2cell type annotation file --------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191004.v1.tsv", data.table = F)

# get the list of cell type marker genes differentially expressed --------------------
renal.markers$gene %>% unique() %>% length()
## 77 genes
genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- c("CA9", genes2plot)

# calculate percent of cells expressing marker genes for each cluster --------------
DefaultAssay(renal.int.nonimmune.obj) <- "RNA"
p <- DotPlot(object = renal.int.nonimmune.obj, features = genes2plot, col.min = 0)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black"),
           strip.text.x = element_text(angle = 90, vjust = 0.5),
           strip.placement = "outside")
p
file2write <- paste0(dir_out, "Dotplot_NonImmune_Marker_Exp.", run_id, ".pdf")
pdf(file2write, width = 16, height = 6)
print(p)
dev.off()



