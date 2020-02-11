# Yige Wu @WashU Dec 2019
## for plotting the HIF target gene expression in pathways

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
# BiocManager::install("pathview")
library(pathview)
packages = c(
  "clusterProfiler",
  "pathview",
  "org.Hs.eg.db"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input the HIF downstream to plot ----------------------------------------
HIF_targets2plot <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/proteomics/plot_HIF_pahtway_protein_heatmap/ccRCC_snRNA_Downstream_Processing - HIF_Target_Summary.tsv", data.table = F)
genes2plot <- HIF_targets2plot$Gene_Symbol[HIF_targets2plot$Tumor == 1 & is.na(HIF_targets2plot$Stromal) & is.na(HIF_targets2plot$Immune)]


filename=system.file("extdata/gse16873.demo", package = "pathview")
filename
gse16873=read.delim(filename, row.names=1)
gse16873.d=gse16873[,2*(1:6)]-gse16873[,2*(1:6)-1]
gse16873.d <- as.matrix(gse16873.d)

gse16873.d.vec <- gse16873.d[,1]
names(gse16873.d.vec) <- rownames(gse16873.d)
data(demo.paths)

data(paths.hsa)

gse16873_d3 <- data.frame(gse16873_d3)

test.gene.data <- as.vector(gse16873_d3[,2])
test.gene.data
names(test.gene.data) <- as.character(gse16873_d3$GeneID)


i <- 1

pv.out <- pathview(gene.data = test.gene.data, pathway.id = "00640",
                   species = "hsa", out.suffix = "gse00640",  gene.idtype = "entrez")

str(pv.out)
head(pv.out$plot.data.gene)

pv.out
