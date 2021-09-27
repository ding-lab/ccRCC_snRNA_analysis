library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(reshape2)
library(ggfortify)
library(tidyverse)



score=read.table('out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
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

annot=score
annot=annot[!duplicated(annot$cell_type),]
rownames(annot)=annot$cell_type
annot=annot[rownames(mydata.cor),]

col_cell_t=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
names(col_cell_t)=c("Tumor","PT","Distal convoluted tubule","Intercalated cells","Loop of Henle","EMT tumor cells","Principle cells","Podocytes")
row_ha=rowAnnotation(Cell_type=annot$cell_type1,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t))


h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F,show_column_dend=FALSE,show_row_dend=FALSE)
pdf('Heatmap_Correlation.NAT_cellTypes.PTbyCluster.2021-08-05.pdf',width=12,height=10)
h
dev.off()
