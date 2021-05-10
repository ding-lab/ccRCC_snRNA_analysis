

# source and libarary -----------------------------------------------------
system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


# input dependencies ------------------------------------------------------
## input the ATAC object
ATAC=readRDS('../15_snATAC_Merged_ccRCC.chromVar.20210219.rds')
###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../Updaate_annotataion/16_snATAC_ccRCC_samples.20210321.tsv',sep='\t',header=TRUE)
## input the motifs
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
                  header=TRUE)

# preprocess --------------------------------------------------------------
## process cell type annotation
cell_t$individual_barcode=rownames(cell_t)
## process ATAC meta data
orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
## reorder the cell type annotation rows
cell_t=cell_t[orig_1$individual_barcode,]
## add Piece_ID, cell type manual
ATAC$Piece_ID=cell_t$Piece_ID
ATAC$cell_type_manual_5=cell_t$cell_type_manual_5
## create test group - divide the tumor cells from different samples
ATAC$test=as.character(ATAC$cell_type_manual_5)
ATAC$test=ifelse(ATAC$test=="Tumor cells",paste(as.character(ATAC$Piece_ID),"Tumor",sep="_"), ATAC$test)
## subset the ATAC object to tumor cells and PT cells from NAT only
atac=ATAC
ATAC=subset(atac, (cell_type_manual_5 %in% c("PT") & Piece_ID %in% c("C3N-01200-N", "C3L-00088-N")) | (cell_type_manual_5=="Tumor cells" & Piece_ID %in% c("C3L-00088-T1", 
                                                                                                                                            "C3L-00917-T1", "C3L-00610-T1", "C3L-01287-T1", "C3L-01313-T1", "C3L-00416-T2", "C3N-00733-T1", 
                                                                                                                                            "C3N-01200-T1", "C3L-00088-T2", "C3L-00079-T1", "C3L-00448-T1","C3L-00096-T1", "C3N-00242-T1",
                                                                                                                                            "C3N-00317-T1")))
## remake test group
ATAC$test=ifelse(as.character(ATAC$Piece_ID) %in% c("C3N-01200-N","C3L-00088-N"),"PT",
                 paste("Tumor",as.character(ATAC$Piece_ID),sep="_"))
## set default assay and get chromVAR assay
DefaultAssay(ATAC) <- 'chromvar'
chromv= GetAssayData(object = ATAC)
## merge motif data with the chromVAR assay data
mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
### remove the extra motif names
res=res[,1:(ncol(res)-2)]
## ?
cell_types=unique(as.character(ATAC$test))
## ?
ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)
## subset and re-order the chromVar assay columns
res=res[,rownames(ann_col0)]


# process by tumor sample -------------------------------------------------
tumor_samples=paste("Tumor",c("C3L-00088-T1",
                              "C3L-00088-T2","C3L-00917-T1","C3L-00448-T1","C3L-00610-T1","C3N-00733-T1","C3L-00079-T1","C3L-01287-T1",
                              "C3N-01200-T1","C3L-00416-T2","C3L-01313-T1","C3L-00096-T1", "C3N-00242-T1","C3N-00317-T1"),sep="_")
final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

for (cell_t1 in tumor_samples){
    print (cell_t1)
    cell_t2='PT'
    ## get the chromVAr assay data for the PT cells
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    ## get the chromVAr assay data for the tumor cells for this tumor sample
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    ## test for each motif
    all_wilcoxon_stat=NULL
    for (motif in 1:nrow(res)){
        ## extract motif scores for each group for this motif
        scores_group1 <- as.numeric(as.character(unlist(res_1[motif,])))
        scores_group2 <- as.numeric(as.character(unlist(res_2[motif,])))
        ## do wilcox test between 2 groups of motif scores
        w_test=wilcox.test(scores_group1, scores_group2)
        ## calculate the mean score of each group
        mean_score1=mean(scores_group1,na.rm=TRUE)
        mean_score2=mean(scores_group2,na.rm=TRUE)
        ## combine results and store
        stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
        all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
    }
    ## reformat the result data frame
    all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
    all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
    all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
    all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
    all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
    all_wilcoxon_stat$diff=abs(all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2)
    all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]
    ### change column names
    colnames(all_wilcoxon_stat)[1:2]=c('cell_t2','TF')
    colnames(all_wilcoxon_stat)[5]=c('pvalue')
    all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
    ## store results
    final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}
## reformat the result data frame
final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)	
colnames(final_wilcoxon_stat)[1:2]=c('cell_t2','TF_Name')
final_wilcoxon_stat$diff=final_wilcoxon_stat$mean_score2-final_wilcoxon_stat$mean_score1
final_wilcoxon_stat=final_wilcoxon_stat[order(-final_wilcoxon_stat$diff),]

# write output -------------------------------------------------------------
write.table(final_wilcoxon_stat,paste("out/Score_difference.Tumor_Normal_comparison.20210319.tsv",sep=""),
            quote=FALSE,sep="\t",row.names=FALSE)
