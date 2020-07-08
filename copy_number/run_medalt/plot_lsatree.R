# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(igraph)
source(paste("./Resources/snRNA_Processed_Data/MEDALT/scripts/github/LSAfunction.R",sep=""))
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ----------------------------------------------------
## input cell tree file
celltree <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/MEDALT/outputs/example_1500cells/CNV.tree.txt")
## input geneLSA
geneLSA <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/MEDALT/outputs/example_1500cells/gene.LSA.txt")
## input bandLSA
bandLSA <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/MEDALT/outputs/example_1500cells/segmental.LSA.txt")

# plot --------------------------------------------------------------------
allsig <- rbind(geneLSA, bandLSA)
LSAnetwork=CNAconnect(allsig,celltree)
nodes=data.frame(id=union(LSAnetwork[,1],LSAnetwork[,2]),size=5)
tab=table(as.character(allsig$cell))
index=match(nodes$id,names(tab))
nodes$size[!is.na(index)]=nodes$size[!is.na(index)]*tab[index[!is.na(index)]]/5
nodes$size[nodes$size<=5]=5
nodes$color="gray"
nodes$color[!is.na(index)]=rainbow(length(unique(allsig$cell)))
annotation=c()
for (i in 1:dim(nodes)[1]){
  if (as.character(nodes$id[i]) %in% as.character(allsig$cell)){
    CNA=allsig[as.character(allsig$cell)==as.character(nodes$id[i]),]
    #pvalue=apply(CNA[,4:5],1,min)
    CNA=CNA[order(CNA$pvalue),]
    CNA=paste(as.character(CNA$region),as.character(CNA$CNA),sep=":")
    CNA1=CNA[1]
    if (length(CNA)>1){
      for (j in 2:min(3,length(CNA))){
        CNA1=paste(CNA1,CNA[j],sep=";")
      }
    }
    annotation[i]=CNA1
  }
}
nodes$annotation=annotation
nodes$size=nodes$size/max(nodes$size)*30
links=data.frame(from=LSAnetwork[,1],to=LSAnetwork[,2],weight=as.numeric(LSAnetwork[,3]))
pdf(file=paste(dir_out,"LSA.tree.pdf",sep=""),width = 10,height = 10,useDingbats = F)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
plot(net, layout=layout_as_tree,vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label.cex=0.5,vertex.label=nodes$annotation)
dev.off()
