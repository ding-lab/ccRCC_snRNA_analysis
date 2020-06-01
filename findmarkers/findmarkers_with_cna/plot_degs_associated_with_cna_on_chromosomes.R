# Yige Wu @WashU May 2020
## run DEG analysis for cells with certain CNV vs CN neutral

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/findallmarkers_genelevel_expected_cnv_vs_neutral_in_tumorcells/20200518.v1/ExpectedCNV_vs_Neutral..FindAllMarkers.Wilcox..20200518.v1.tsv")
## get genes in the txDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# transform gene symbol to entrez idss ----------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## filter markers
markers_df <- markers_df %>%
  filter(p_val_adj < 0.05)
## retrieve UCSC stable ids
genesymbol2ucsc_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                            filters = 'hgnc_symbol', 
                            values = unique(markers_df$de_gene_symbol), 
                            mart = ensembl)
## add entrez ids to the deg table
markers_df$entrezgene_id <- mapvalues(x = markers_df$de_gene_symbol, from = genesymbol2ucsc_df$hgnc_symbol, to = as.vector(genesymbol2ucsc_df$entrezgene_id))
## filter markers by entrez id
markers_df <- markers_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != de_gene_symbol)

# process by cna gene -----------------------------------------------------
for (cna_gene_plot in unique(markers_df$cna_gene_symbol)) {
  ## set cna expected state
  # cna_gene_plot <- "VHL"
  cna_gene_expected_state <- unique(markers_df$cna_gene_state.1[markers_df$cna_gene_symbol == cna_gene_plot])
  
  # get gene position from txdb ---------------------------------------------
  hs.genes <- genes(txdb, columns = c("gene_id"))
  head(hs.genes)
  names(hs.genes)
  
  # map genes to genome position --------------------------------------------
  ## get genes
  markers_for_cnagene_df <- markers_df %>%
    filter(cna_gene_symbol == cna_gene_plot) %>%
    group_by(de_gene_symbol, entrezgene_id) %>%
    summarise(num_aliquot_de = n(), mean_avg_logFC = mean(avg_logFC))
  markers_for_cnagene_df <- as.data.frame(markers_for_cnagene_df)
  ## make entrez id the row names
  rownames(markers_for_cnagene_df) <- markers_for_cnagene_df$entrezgene_id
  
  ## add differential expression result to the GRange object
  mcols(hs.genes) <- markers_for_cnagene_df[names(hs.genes), c("mean_avg_logFC", "entrezgene_id", "de_gene_symbol", "num_aliquot_de")]
  head(hs.genes, n=4)
  ## filter GRange object
  filtered.hs.genes <- hs.genes[!is.na(hs.genes$mean_avg_logFC)]
  head(filtered.hs.genes, n=4)
  ## order GRange object
  ordered <- filtered.hs.genes[order(filtered.hs.genes$num_aliquot_de, na.last = TRUE, decreasing = T),]
  head(ordered, n=20)
  ## get top genes
  num_top_genes <- length(which(ordered$num_aliquot_de == max(ordered$num_aliquot_de)))
  num_top_genes <- max(c(num_top_genes, 20))
  num_top_genes <- min(c(num_top_genes, length(ordered)))
  top.genes <- ordered[1:num_top_genes]
  top.genes
  # plot --------------------------------------------------------------------
  ## get plot parameters
  fc.ymax <- signif(max(abs(range(filtered.hs.genes$mean_avg_logFC))), digits = 2)
  fc.ymin <- -fc.ymax
  cex.val <- filtered.hs.genes$num_aliquot_de/2
  points.top <- 0.8
  ### make color
  col.over <- "#FFBD07AA"
  col.under <- "#00A6EDAA"
  sign.col <- rep(col.over, length(filtered.hs.genes))
  sign.col[filtered.hs.genes$mean_avg_logFC<0] <- col.under
  ### parameters for ideogram
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 10
  ## plot
  file2write <- paste0(dir_out, cna_gene_plot, "_", cna_gene_expected_state, "_associated_degs.pdf")
  pdf(file = file2write, width = 8, height = 15)
  kp <- plotKaryotype(genome = "hg38")
  kpAddMainTitle(kp, main = paste0("Differentially Expressed Genes in Cells with ", cna_gene_plot, " ", cna_gene_expected_state, " vs Neutral"))
  kpPoints(kp, data=filtered.hs.genes, y=filtered.hs.genes$mean_avg_logFC, ymax=fc.ymax, ymin=fc.ymin, cex = cex.val, r1=points.top, col=sign.col)
  kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
  gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
  kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$mean_avg_logFC, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=0.5, col="#777777")
  kpPlotMarkers(kp, top.genes, labels = top.genes$de_gene_symbol, text.orientation = "horizontal", r0=0.5, label.dist = 0.008, label.color="#444444", line.color = "#777777")
  # kp <- kpPlotDensity(kp, data=filtered.hs.genes, window.size = 10e4, data.panel = 2)
  dev.off()
}



