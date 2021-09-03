# Yige Wu @WashU Aug 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
score=read.table('./Resources/snATAC_Processed_Data/Enriched_Motifs//Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)

# process data matrix-----------------------------------------------------------------
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1=score_1[!(score_1$cell_type1 %in% c(12,14)),]
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:14),"PT",score_1$cell_type1)
score_1=score_1[!is.na(score_1$cell_type1),]

score=score_1
score=score[!(score$cell_type1 %in% c("Macrophages","NK cells","Unknown","CD8+ T-cells","CD4+ T-cells","DC","Endothelial cells","Fibroblasts")),]
score$cell_type1=ifelse(score$cell_type1=="Microglia","TAM",score$cell_type1)
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

mydata.cor = cor(t(all_2), method = c("spearman"))

# make colors -------------------------------------------------------------
## colors for the cell types
colors_celltype1 <- colors_cellgroup14[c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts")]
names(colors_celltype1)[1] <- "Tumor"
colors_celltype2 <- c(colors_cellgroup14[c("EMT tumor cells", "Normal epithelial cells", "Myofibroblasts", "Immune others", "B-cells")], 
                 Polychrome::palette36.colors(n = 36)[c("Vivid_Violet","Light_Olive_Brown", "Very_Light_Blue")], "grey80")
names(colors_celltype2) <- c("EMT tumor cells", "PT", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', 
                        "Intercalated cells", "Podocytes", "Endothelial cells", "Unknown")
colors_celltype <- c(colors_celltype1, colors_celltype2)
colors_celltype["Distal convoluted tubule"] <- brewer.pal(n = 8, name = "Paired")[8]
colors_celltype["EMT tumor cells"] <- brewer.pal(n = 8, name = "Paired")[6]
colors_celltype <- colors_celltype[unique(annot$cell_type1)]
## colors for the heatmap
colors_heatmapbody <- colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))

# make row annotation -----------------------------------------------------
annot=score
annot=annot[!duplicated(annot$cell_type),]
rownames(annot)=annot$cell_type
annot=annot[rownames(mydata.cor),]
annot <- annot %>%
  mutate(subcluster = ifelse(cell_type1 == "PT", paste0("PT_C", Piece_ID), ""))
left_row_anno=rowAnnotation(subcluster = anno_mark(at = which(annot$cell_type1 == "PT"), 
                                                   labels = annot$subcluster[annot$cell_type1 == "PT"], side = "left",
                                                   labels_gp = gpar(fontsize = 25)),
                            Cell_type=anno_simple(annot$cell_type1, col = colors_celltype, width = unit(0.75, "cm")), show_legend = F)

# plot --------------------------------------------------------------------
h=Heatmap(matrix = mydata.cor, 
          col= colors_heatmapbody, show_heatmap_legend = F,
          ## column
          show_column_names = FALSE, show_column_dend=FALSE,
          ## row 
          left_annotation=left_row_anno, 
          # right_annotation = right_row_anno,
          show_row_names=F,show_row_dend=FALSE)
list_lgd = list(
  Legend(title = "Correlation of TF\nmotif scores", title_gp = gpar(fontsize = 20),
         col_fun = colors_heatmapbody, labels_gp = gpar(fontsize = 20), legend_height = unit(4, "cm"), grid_width = unit(0.75, "cm")),
  Legend(labels = names(colors_celltype), labels_gp = gpar(fontsize = 20),
         title = "Cell type", title_gp = gpar(fontsize = 20), grid_width = unit(0.75, "cm"),
         legend_gp = gpar(fill = colors_celltype), border = NA))


file2write <- paste0(dir_out, 'Heatmap_Correlation.NAT_cellTypes.PTbyCluster.2021-08-05.pdf')
pdf(file2write,width=14,height=10, useDingbats = F)
draw(object = h,
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()



