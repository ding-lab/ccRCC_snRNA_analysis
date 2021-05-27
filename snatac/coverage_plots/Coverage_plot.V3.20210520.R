#system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")
###For some samples (with many cells >6K) python-package used in chromVar doesn't work properly; Need to use this 2 commands:
###export OMP_NUM_THREADS=1
###export USE_SIMPLE_THREADED_LEVEL3=1
###export OPENBLAS_NUM_THREADS=1

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
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyverse)

ATAC=readRDS(paste('../3.Merge_snATAC/Merge.SelectPeaks.v.20210503/',
'26_ccRCC_snATAC.selectedPeaks.chromvar.CICERo.v6.20210512.rds',sep=''))

Idents(ATAC)=ATAC$Piece_ID


hif1_a_genes=c('PGK1', 'LDHA', 'NOL3', 'EGLN3', 'FAM162A', 'BNIP3', 'HK1', 'VEGFA', 'TGFA')
other_gene=c('CA9', 'VIM', 'EGFR', 'PKM', 'PFKP', 'GAPDH', 'ENO2', 'TPI1', 'SLC2A1', 'BHLHE41')


atac=ATAC
ATAC=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% c('C3L-00004-T1', 'C3L-00010-T1',
 'C3L-00079-T1', 'C3L-00088-T2', 'C3L-00088-T1',  'C3L-00448-T1',
'C3L-00583-T1', 'C3L-00610-T1', 'C3L-00790-T1', 'C3L-00908-T1', 'C3L-00917-T1',
'C3L-01287-T1',  'C3L-01313-T1', 'C3N-00317-T1','C3N-00733-T1',
'C3N-01200-T1','C3N-01213-T1','C3L-00416-T2','C3L-01302-T1','C3L-00096-T1', 'C3L-00026-T1', 
'C3N-00242-T1', 'C3N-00437-T1', 'C3N-00495-T1')) | cell_type=='PT' & Piece_ID %in%
c('C3L-00088-N','C3N-01200-N'))

Idents(ATAC)=factor(ATAC$Piece_ID,levels=c('C3L-00004-T1', 'C3L-00010-T1',
 'C3L-00079-T1', 'C3L-00088-T2', 'C3L-00088-T1',  'C3L-00448-T1',
'C3L-00583-T1', 'C3L-00610-T1', 'C3L-00790-T1', 'C3L-00908-T1', 'C3L-00917-T1',
'C3L-01287-T1',  'C3L-01313-T1', 'C3N-00317-T1','C3N-00733-T1',
'C3N-01200-T1','C3N-01213-T1','C3L-00416-T2','C3L-01302-T1','C3L-00096-T1', 'C3L-00026-T1', 
'C3N-00242-T1', 'C3N-00437-T1', 'C3N-00495-T1', 
'C3L-00088-N','C3N-01200-N'))

Idents(ATAC)=factor(ATAC$cell_type,levels=c('Tumor','PT'))

#######
peaks_list=c('chr6-44082757-44083257', 'chr6-43768858-43769358')

#degs=read.table('DEG_associated_Peaks.Promoter.Prioritized.20210518.v1.tsv',sep='\t',header=T)
peaks=read_delim(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'Signac.1.0.0/6.DA_motifs/Match_motifs/Motifs_matched.DEG_associated_Peaks.20210517.v1.tsv',
sep=''),delim='\t')
peaks=as.data.frame(peaks)
peaks_s=peaks[peaks$Peak %in% peaks_list,]
#peaks_s=peaks[peaks$Peak_Type=='Promoter' & peaks$Type=='Promoter' & peaks$Gene==peaks$Gene_1,]


#peaks_s=peaks_s[peaks_s$Gene %in% c(degs$Gene),]

peaks=peaks_s
peaks=peaks[!duplicated(peaks[c('Peak','TF_name')]),]

#VEGFA_peaks.20210520

Idents(ATAC)=factor(ATAC$cell_type,levels=c('Tumor','PT'))

unique_p=unique(peaks$Peak)
for (i in 1:length(unique_p)){
peak=unique_p[i]
chr=gsub('(.*)-(.*)-(.*)','\\1',peak)
st=as.numeric(gsub('(.*)-(.*)-(.*)','\\2',peak))
en=as.numeric(gsub('(.*)-(.*)-(.*)','\\3',peak))
new_st=st-1000
new_en=en+1000
new_peak=paste(chr,new_st,new_en,sep='-')

peaks_g=peaks[peaks$Peak==peak,]
tfs=c('ARNT::HIF1A','EGR1','HIF1A','HSF2','KLF14','KLF15','MXI1','NEUROD1','NEUROG2(var.2)',
'NFKB1','NFKB2','NRF1','RBPJ','SP1','SP3','SP9','SREBF2','TCFL5','ZBTB14','ZNF148','ZNF75D')
peaks_g=peaks_g[peaks_g$TF_name %in% tfs,]
gene=unique(peaks_g$Gene)
for (tf in peaks_g$TF_name){
peaks_tf=peaks_g[peaks_g$TF_name ==tf,]
if (nrow(peaks_tf)>0){
p=CoveragePlot(
  object = ATAC,
  region = new_peak,
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(peaks_tf$motif_coords, sep = c("-", "-")),
  ranges.title = paste(tf,"-motifs",sep=''),
  links=FALSE
)
dir.create(paste('VEGFA_peaks.20210520/',tf,sep=''))
pdf(paste("VEGFA_peaks.20210520/",tf,"/",gene,"_",tf,"_",new_peak,".pdf",
sep=""),width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()
}
print(paste(new_peak,i,gene,sep=' '))
}
}

##############################
#####plot by Sample
###############


unique_p=unique(peaks$Peak)
for (i in 1:length(unique_p)){
peak=unique_p[i]
chr=gsub('(.*)-(.*)-(.*)','\\1',peak)
st=as.numeric(gsub('(.*)-(.*)-(.*)','\\2',peak))
en=as.numeric(gsub('(.*)-(.*)-(.*)','\\3',peak))
new_st=st-1000
new_en=en+1000
new_peak=paste(chr,new_st,new_en,sep='-')

peaks_g=peaks[peaks$Peak==peak,]
tfs=c('ARNT::HIF1A','EGR1','HIF1A','HSF2','KLF14','KLF15','MXI1','NEUROD1','NEUROG2(var.2)',
'NFKB1','NFKB2','NRF1','RBPJ','SP1','SP3','SP9','SREBF2','TCFL5','ZBTB14','ZNF148','ZNF75D')
peaks_g=peaks_g[peaks_g$TF_name %in% tfs,]
gene=unique(peaks_g$Gene)
for (tf in peaks_g$TF_name){
peaks_tf=peaks_g[peaks_g$TF_name ==tf,]
if (nrow(peaks_tf)>0){
p=CoveragePlot(
  object = ATAC,
  region = new_peak,
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(peaks_tf$motif_coords, sep = c("-", "-")),
  ranges.title = paste(tf,"-motifs",sep=''),
  links=FALSE
)
dir.create(paste('VEGFA_peaks.BySample.20210520/',tf,sep=''))
pdf(paste("VEGFA_peaks.BySample.20210520/",tf,"/",gene,"_",tf,"_",new_peak,".pdf",
sep=""),width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()
}
print(paste(new_peak,i,gene,sep=' '))
}
}
#########


####Trye to make CICERO PLOT:
####Links are absent?
conns=read.table('../7.CICERO/out/26_ccRCC_snATAC_CICERO.tsv',sep='\t',header=TRUE)
ccans=read.table('../7.CICERO/out/26_ccRCC_snATAC_CICERO.CCAN.tsv',sep='\t',header=TRUE)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(ATAC) <- links

####To many links... try different way:
#conns_1=conns[conns$Peak1 %in% unique_p | conns$Peak2 %in%]
conns_1=conns[conns$coaccess>0.5,]
conns_1=conns_1[!is.na(conns_1$Peak1),]
ccans_1=ccans[ccans$Peak %in% c(conns_1$Peak1,conns_1$Peak2),]
links_1<- ConnectionsToLinks(conns = conns_1, ccans = ccans_1)
Links(ATAC)<-links_1


p=CoveragePlot(
  object = ATAC,
  region ='chr6-43766858-44085257',
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(unique_p, sep = c("-", "-")),
  ranges.title = paste("VEGFA_peaks",sep=''),
  links=TRUE
)
pdf('VEGFA_peaks.BySample.20210520/Links_VEGFA.pdf',width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()
pdf('VEGFA_peaks.20210520/Links_Tumor_PT_VEGFA.pdf',width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()










add=PeakPlot(
  ATAC,
  region=new_peak,
  color = "dimgrey",
  peaks=StringToGRanges(peaks_tf$motif_coordsc, sep = c("-", "-"))
)

x=CombineTracks(
  plotlist = list(p, add))













########################
###V2, without 1287:####
########################

peaks=read.table('BAP1_DOWN_top.20210408.witout1287.tsv',sep='\t',header=TRUE)
peaks_1=StringToGRanges(peaks$peak, sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$peak_distanceToTSS=anno$distanceToTSS
peaks=peaks[peaks$Type=="Promoter",]
chrs=gsub('(.*)-.*-.*','\\1',peaks$peak)
table(chrs)

for (i in 1:nrow(peaks)){
peak=peaks$peak[i]
chr=gsub('(.*)-(.*)-(.*)','\\1',peak)
st=as.numeric(gsub('(.*)-(.*)-(.*)','\\2',peak))
en=as.numeric(gsub('(.*)-(.*)-(.*)','\\3',peak))
new_st=st-1000
new_en=en+1000
new_peak=paste(chr,new_st,new_en,sep='-')
gene=peaks$Gene[i]
p=CoveragePlot(
  object = ATAC,
  region = new_peak,
  annotation = TRUE,
  peaks = TRUE
)


pdf(paste("coverage_plots_without1287/",gene,"_",new_peak,".pdf",sep=""),width=8,height=8,useDingbats=FALSE)
print(p)
dev.off()
print(paste(new_peak,i,gene,sep=' '))
}







CoveragePlot(
  object = ATAC,
  region = "chr11-119094215-119096715",
  annotation = TRUE,
  peaks = TRUE,
)
CoveragePlot(
  object = ATAC,
  region = "chr14-103714213-103716713",
  annotation = TRUE,
  peaks = TRUE,
)


####The plot looks as expected:
res_1=res
res_1=res_1[,c(1:5,7:8)]
colnames(res_1)[1:5]=c('Case','Cell_count','Total_fragments','MGMT_promoter_fragm','MGMT_NormalizedByTotal')
write.table(res_1, "MGMT_status_byATAC_accessibility.tsv",sep='\t',row.names=FALSE,quote=FALSE)








###################################################################
###Previousattmepts,worse correspondence with the coverage plot:###
###################################################################

y1=y[y$Norm_cov!=0,]
st=aggregate(y1$Norm_cov, by=list(y1$Case),FUN='mean')
st=st[order(st$x),]


subt=read.table(paste(data_dir,'akshay_p300_acetyl_status.tsv',sep=''),sep='\t',header=TRUE)

meta=read.table(paste(data_dir,'13_samples_Annotated.v.20210201.ManuallyReiewed.tsv',sep=''),
sep='\t',header=TRUE)

orig_1=as.data.frame(ATAC$dataset)
orig_1$i_barc=rownames(orig_1)

meta=meta[orig_1$i_barc,]

ATAC$cell_type_manual=meta$cell_type_manual
ATAC$Piece_ID=meta$Piece_ID

#atac=ATAC
ATAC=subset(ATAC, cell_type_manual=="Tumor" & Piece_ID %in% c('C3L-02705', 'C3N-01334', 'C3N-01518', 'C3N-01816',
 'C3N-01818', 'C3N-02186', 'C3N-02188', 'C3N-02769', 'C3N-02784', 'C3N-03186', 'C3N-00663'))

mgmt=read.table(paste(data_dir,"mgmt_cpg_island_dna_methyl_probes.tsv",sep=''),sep='\t',header=TRUE)
bady_probes=c('cg12434587','cg12981137')
wenwei_probes = c('cg00618725', 'cg12434587', 'cg12981137')
more_bady_2012_probes = c('cg14194875', 'cg00618725', 'cg12434587', 'cg16215402', 'cg18026026',
'cg12981137')

mgmt$regions=paste(mgmt$seqnames,mgmt$start,mgmt$end,sep='-')
p_bady=mgmt[mgmt$probe_id %in% bady_probes,]
p_ww=mgmt[mgmt$probe_id %in% wenwei_probes,]
p_more=mgmt[mgmt$probe_id %in% more_bady_2012_probes,]

probes_bady=StringToGRanges(regions = p_bady$regions)
probes_ww=StringToGRanges(regions = p_ww$regions)
probes_more=StringToGRanges(regions = p_more$regions)
probes=StringToGRanges(regions = mgmt$regions)

Idents(ATAC)=ATAC$Piece_ID
#Idents(ATAC)=factor(Idents(ATAC),levels=c('C3L-02705', 'C3N-01518', 'C3N-02769', 'C3N-00663', 'C3N-02784', 
#'C3N-01334', 'C3N-01816', 'C3N-01818', 'C3N-02186', 'C3N-02188', 'C3N-03186'))

Idents(ATAC)=factor(Idents(ATAC),levels=st$Group.1)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(ATAC) <- annotations


CoveragePlot(
  object = ATAC,
  region = "chr10-129465872-129468372",
#  features = "MGMT",
  annotation = TRUE,
  peaks = TRUE,
#  tile = TRUE,
  links = TRUE
)

#DotPlot(object=ATAC,features='MAPK8'

region='chr10-129466872-129467372'
peak_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey"
)
probe_all_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes
)
probe_bady_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_bady
)
probe_ww_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_ww
)
probe_more_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_more
)
gene_plot <- AnnotationPlot(
  object = ATAC,
  region = region
)

cov_plot<-CoveragePlot(
  object = ATAC,
  region = region,
  annotation = FALSE,
  peaks = FALSE
)

CombineTracks(
  plotlist = list(cov_plot, peak_plot, probe_bady_plot,probe_ww_plot,probe_more_plot, gene_plot)
)


####V2:
CombineTracks(
  plotlist = list(cov_plot, peak_plot, probe_bady_plot,probe_ww_plot,probe_more_plot,probe_all_plot, gene_plot)
)





subt=read.table(paste(data_dir,'gbm_all_subtype_collections.v5.1.tsv',sep=''),sep='\t',header=TRUE)
