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

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)


## input differential expression analysis
deg_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)

# input TF table ----------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
HIF_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1"))

# specify the HIF downstream to plot ----------------------------------------
# HIF_targets2plot <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/proteomics/plot_HIF_pahtway_protein_heatmap/ccRCC_snRNA_Downstream_Processing - HIF_Target_Summary.tsv", data.table = F)
# genes2plot <- HIF_targets2plot$Gene_Symbol[HIF_targets2plot$Tumor == 1 & is.na(HIF_targets2plot$Stromal) & is.na(HIF_targets2plot$Immune)]
genes2plot <- HIF_tf_tab$target_genesymbol
genes2plot <- unique(c("PLIN2", "CA12", "PLG", "IL6", "ADM", "VEGFA", "BNIP3", "HK1", "HK2", "PFK", "ALDOA", "PGK1", "LDHA", "NOS2", "ABL2", "EPO", "POUSF1", "SCGB3A1", "TGFA", "CCND1", "DLL4", "ANGPT2",
                       genes2plot))

# convert gene symbol to entrez id ----------------------------------------
entrez_ids = bitr(genes2plot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# KEGG pathway mapping, the only one supported by pathview------------------------------------------------------
kk <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk_result <- kk@result

# combine gene symbol, entrez ids with fold change -------------------------------------
## get the average of the average fold change among tumor cells
fc_tab <- deg_tab %>%
  filter(gene %in% genes2plot) %>%
  filter(cluster %in% cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"]) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(max_avg_logFC = max(avg_logFC, na.rm = T))

## map values
fc_tab$entrez_id <- mapvalues(x = fc_tab$gene, from = entrez_ids$SYMBOL, to = entrez_ids$ENTREZID)

# pathview ----------------------------------------------------------------
geneList <- fc_tab$max_avg_logFC
names(geneList) <- fc_tab$entrez_id
hsa04066 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04066",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
hsa00010 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa00010",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))



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
