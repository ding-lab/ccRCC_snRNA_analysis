# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(igraph)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input cell tree file ----------------------------------------------------
celltree <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/MEDALT/outputs/example_1500cells/CNV.tree.txt")


# Plot cell tree ----------------------------------------------------------
nodes=data.frame(id=union(as.character(celltree[,1]),as.character(celltree[,2])),size=3)
nodes$color="lightblue"
nodes$color[nodes$id==setdiff(as.character(celltree[,1]),as.character(celltree[,2]))]="black"
nodes <- rbind(nodes[nodes$color != "black",], nodes[nodes$color == "black",])
net <- graph_from_data_frame(d=celltree, vertices=nodes, directed=T)
pdf(file=paste(dir_out,"singlecell.tree.pdf",sep=""),width = 20,height = 20,useDingbats = F)
plot(net, vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label=NA)
dev.off()
