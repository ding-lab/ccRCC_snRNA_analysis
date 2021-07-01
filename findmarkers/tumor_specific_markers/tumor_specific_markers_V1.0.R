################################################################################
#######################cell surface markers DE genes ###########################
#######################           Lijun Yao         ############################
#######################          2020-09-03        ###########################
################################################################################
library(Seurat)
library(dplyr)
library(getopt)
library(biomaRt)
library(optparse)

option_list = list(
  make_option(c("-o", "--output_path"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-r","--rds_path"),
              type="character",
              default=NULL,
              help="path of a file with the list of sample names and rds paths",
              metavar="character"),
  make_option(c("-c","--cell_type"),
              type="character",
              default=NULL,
              help="cell type of interest",
              metavar="character"),
  make_option(c("-s","--cspa"),
              type="character",
              default=NULL,
              help="cellular location from Cell Surface Protein Atlas",
              metavar="character"),
  make_option(c("-p","--hpa"),
              type="character",
              default=NULL,
              help="cellular location from Human Protein Atlas",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rds_path)){
  print_help(opt_parser)
  stop("Path to data is required (--rds_path).n", call.=FALSE)
}
if (is.null(opt$output_path)){
  print_help(opt_parser)
  stop("Path to output folder is required (--output_path).n", call.=FALSE)
}
if (is.null(opt$cell_type)){
  print_help(opt_parser)
  stop("cell type of interest is required (--cell_type).n", call.=FALSE)
}
if (is.null(opt$cspa)){
  print_help(opt_parser)
  stop("Path to cellular location file from Cell Surface Protein Atlas is required (--cspa).n", call.=FALSE)
}
if (is.null(opt$hpa)){
  print_help(opt_parser)
  stop("Path to cellular location file from Human Protein Atlas is required (--hpa).n", call.=FALSE)
}

# read in initial arguments
rds_list <- opt$rds_path
out_path <- opt$output_path
tumor_ct <- opt$cell_type
cspa_file <- opt$cspa
hpa_file <- opt$hpa

# make output  make output directory
dir.create(out_path,showWarnings = FALSE)

########################################################################
###step 1: DE analysis function for a given cell type vs other cell types
########################################################################
tumor_vs_rest_DE_fun <- function(sobj,tumor_ct,sample){
  DefaultAssay(sobj)<-"RNA"
  DE_genes <- FindMarkers(sobj,ident.1=tumor_ct) 
  DE_genes$gene_symbol <- rownames(DE_genes)
  DE_genes$sample_id <- rep(sample,nrow(DE_genes))
  rownames(DE_genes) <- paste0(sample_id,".",rownames(DE_genes))
  return(DE_genes)
}

#####################################################################
###step 2: pairwise DE analysis function for compairing between a given cell type and each other cell type
#####################################################################
pairwise_DE_fun <- function(sobj,tumor_ct,sample){
  cell.types=as.character(subset(as.data.frame(table(Idents(sobj))),Freq>=3)$Var1)
  cat("cell types include: ",paste(cell.types,sep=","),"\n")
  non_tumor_ct <- cell.types[!(cell.types %in% c(tumor_ct,"unknown"))]
  pairwise_DE <- as.data.frame(matrix(nrow=1,ncol=8))
  colnames(pairwise_DE) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","sample_id","cell_type","gene_symbol")
  rownames(pairwise_DE) <- "tmp_row"
  DefaultAssay(sobj)<-"RNA"
  for (ct in non_tumor_ct) {
    DE_tmp <- FindMarkers(sobj, ident.1 =tumor_ct, ident.2 = ct)
    DE_tmp$sample_id <- rep(sample,nrow(DE_tmp))
    DE_tmp$cell_type <- rep(ct,nrow(DE_tmp))
    DE_tmp$gene_symbol <- rownames(DE_tmp)
    rownames(DE_tmp) <- paste0(DE_tmp$sample_id,".",DE_tmp$cell_type,".",rownames(DE_tmp))
    pairwise_DE=rbind(pairwise_DE,DE_tmp)
  }
  pairwise_DE <- pairwise_DE[rownames(pairwise_DE)!="tmp_row",]
  return(pairwise_DE)
}


################################################################
###Merging DE genes from all samples
#####################################################################
#generate initial dataframe for merge
tumor_vs_rest_DE_total <- as.data.frame(matrix(nrow=1,ncol=8))
colnames(tumor_vs_rest_DE_total) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","gene_symbol","sample_id","avg_norm_exp")
rownames(tumor_vs_rest_DE_total) <- "tmp_row"
  
pairwise_DE_total<- as.data.frame(matrix(nrow=1,ncol=8))
colnames(pairwise_DE_total) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","sample_id","cell_type","gene_symbol")
rownames(pairwise_DE_total) <- "tmp_row"

# For each sample, read object and run step1, step2 
# rds_list <- "/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/rds_list.txt"
rds_list_df <- read.table(rds_list,header=FALSE,stringsAsFactors=FALSE,sep="\t")
samples <- rds_list_df$V1

#samples <- list.dirs(path ="/diskmnt/Datasets/mmy_scratch/lyao/MMY/Organized/",full.names=FALSE, recursive = FALSE)
#samples <- samples[!(samples %in% c("ND_083017","ND_090617","Normal_sorted_170531","Normal_sorted_170607"))]
#no "25183","57075_2" no tumor cells/all tumor cells, "MMY67868" some clusers not assigned
#tumor_ct="Plasma"
for (sample_id in samples) {
  cat(paste0("READING SEURAT OBJECT FOR SAMPLE ",sample_id,"...\n"))
  sobj <- readRDS(file = as.vector(subset(rds_list_df,V1==sample_id)$V2))
  #skip the sample if it only has the tumor cells or no any tumor cells
  cell.types = as.character(subset(as.data.frame(table(Idents(sobj))),Freq>=3)$Var1)
  if ((length(cell.types)==1) & (tumor_ct %in% cell.types)) {
    cat(paste0(sample," Only Contain Tumor Cells So Skip This Sample"))
    next
  }
  if (!(tumor_ct %in% cell.types)) {
    cat(paste0(sample," Does Not Have Any Tumor Cells So Skip This Sample"))
    next
  }
  
  #step1 tumore vs all the rest cell types
  cat("STEP1 DEG analysis between tumor cells and other cell populations as as whole...\n")
  DE_genes_tmp <- tumor_vs_rest_DE_fun(sobj=sobj,tumor_ct=tumor_ct,sample = sample_id)
  #append avg exp of tumor cells to the DE_genes_tmp df
  row_genes <- rownames(DE_genes_tmp) %>% strsplit("[.]") %>% lapply("[",2) %>% unlist 
  DefaultAssay(sobj)<-"RNA"
  genes <- intersect(row_genes,GetAssayData(object = sobj) %>% rownames)
  norm_exp <- GetAssayData(object = sobj)[genes,WhichCells(sobj, idents = tumor_ct)]
  norm_exp_avg <- apply(norm_exp,1,mean) %>% as.data.frame
  colnames(norm_exp_avg) <- "avg_norm_exp"
  norm_exp_avg$gene_symbol <- rownames(norm_exp_avg)
 
  DE_genes_tmp <- merge(DE_genes_tmp,norm_exp_avg,by="gene_symbol",all.x=TRUE,incomparables = NA,sort=FALSE)
  rownames(DE_genes_tmp) <- paste0(DE_genes_tmp$sample_id,".",DE_genes_tmp$gene_symbol)
  tumor_vs_rest_DE_total <- rbind(tumor_vs_rest_DE_total,DE_genes_tmp)
   
  #step2 tumore vs each other cell type -> pairwise comparison
  cat("STEP2 DEG analysis between tumor cells and each other cell population...\n")
  pairwise_DE_tmp <- pairwise_DE_fun(sobj=sobj,tumor_ct=tumor_ct,sample=sample_id)
  pairwise_DE_total <- rbind(pairwise_DE_total,pairwise_DE_tmp)
}

tumor_vs_rest_DE_total <- tumor_vs_rest_DE_total[rownames(tumor_vs_rest_DE_total)!="tmp_row",]
pairwise_DE_total <- pairwise_DE_total[rownames(pairwise_DE_total)!="tmp_row",]

pairwise_DE_total$cell_type <- sub(" ", "_", pairwise_DE_total$cell_type)
rownames(pairwise_DE_total) <- sub(" ", "_",rownames(pairwise_DE_total))

cat("SAVING THE DEG ANALYSIS RESULTS TO TEXT FILES...\n")
write.table(tumor_vs_rest_DE_total,paste0(out_path,"/",tumor_ct,"_vs_combined_others_DE.txt"),sep="\t",col.names=TRUE,row.names=TRUE,quote = FALSE)
write.table(pairwise_DE_total,paste0(out_path,"/",tumor_ct,"_vs_others_pairwise_DE.txt"),sep="\t",col.names=TRUE,row.names=TRUE,quote = FALSE)

##########################################################################
###step3: filtering
##########################################################################
p_val_adj_filter <- 0.05
avg_logFC_filter <- 0
################################
#########EXAMPLE################
###             count FC>0 p<0.05?
###PMEL 47499   5     4    3       <- will be removed
###PMEL 47491   3     3    3       <- will be kept 
#################################
#lower the threshold, allow one more cell types to be significant highly expressed with plasma
tmp <- pairwise_DE_total %>% dplyr::select(sample_id,gene_symbol,avg_logFC,p_val_adj) %>% dplyr::group_by(sample_id,gene_symbol) %>% dplyr::summarise(significance_count=sum(p_val_adj<p_val_adj_filter),pos_logFC_count=sum(avg_logFC>0),count=length(sample_id)) %>% filter(significance_count==count & pos_logFC_count==count) %>% as.data.frame
tmp$gene_symbol=factor(tmp$gene_symbol,levels=unique(tmp$gene_symbol))

gene_threshold <- unname(quantile(table(tmp$gene_symbol),0.50)) # for a certain marker, it must be reported in more than 50% of samples
DE_genes_filter_1 <- names(table(tmp$gene_symbol))[table(tmp$gene_symbol)>=gene_threshold]

###FILTER 2: Compare tumor cells with other cell types as a whole for each sample
## filter DE_genes_filter_1 by the presence within DE gene list(when compare plasma vs non plasma)
tmp <- tumor_vs_rest_DE_total %>% dplyr::select(gene_symbol,avg_logFC,p_val_adj) %>% dplyr::group_by(gene_symbol) %>% dplyr::summarise(significance_count=sum(p_val_adj<p_val_adj_filter),pos_logFC_count=sum(avg_logFC>avg_logFC_filter),count=length(gene_symbol)) %>% as.data.frame
tmp <- tmp %>% dplyr::filter(pos_logFC_count/count> 0.90 & significance_count/count> 0.75) 
sample_count_threshold <- unname(quantile(tmp$count,0.75)) #threshold top 25%
DE_genes_filter_2 <- tmp %>% dplyr::filter(count>=sample_count_threshold) %>% .$gene_symbol %>% as.character

### taking the intersection of the above two filtering
DE_genes_filtered <- intersect(DE_genes_filter_1,DE_genes_filter_2)
cat(paste0(length(DE_genes_filtered)," genes were found to be tumor specific\n"))

### sort genes by avg_logFC: summary table is from DEG between tumor cells and other cell populations as as whole
tmp <- tumor_vs_rest_DE_total %>% dplyr::filter(gene_symbol %in% DE_genes_filtered) %>% dplyr::group_by(gene_symbol) %>% dplyr::summarise(avg_logFC=mean(avg_logFC)) %>% as.data.frame
DE_genes_filtered <- tmp[rev(order(tmp$avg_logFC)),] %>% .$gene_symbol 
write.table(DE_genes_filtered,paste0(out_path,"/",tumor_ct,"_specific_DEG.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote = FALSE)
#DE_genes_filtered_highranking<-tmp[rev(order(tmp$avg_logFC)),] %>% filter(avg_logFC>0.75) %>% .$gene_symbol %>% as.character

######################################################################
###cell surface markers(Annotate genes with GO term plasma memberane)
#######################################################################
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
plasma_membrane_genes <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), filters = 'go', values = 'GO:0005886', mart = ensembl)$hgnc_symbol
DE_genes_filtered_surface <- intersect(DE_genes_filtered,plasma_membrane_genes)
cat(paste0(length(DE_genes_filtered_surface)," genes were found on the cell surface based on GO term\n"))

#####################################################################################################################
###Append the information of FC and adjusted p value to the DE_genes_filtered data frame, p_val_adj and avg_logFC are from DE analysis of tumor vs the rest as a whole
#####################################################################################################################
DE_genes_filtered_df <- as.data.frame(DE_genes_filtered)

fc_vector <- c()
fdr_vector <- c()
exp_vector <- c()
freq_vector <- c()
sample_vector <- c()
for (gene in DE_genes_filtered_df$DE_genes_filtered){
  fc <- mean(subset(tumor_vs_rest_DE_total,gene_symbol==gene)$avg_logFC)
  fdr <- mean(subset(tumor_vs_rest_DE_total,gene_symbol==gene)$p_val_adj)
  exp <- mean(subset(tumor_vs_rest_DE_total,gene_symbol==gene)$avg_norm_exp)  
  freq <- length(unique(subset(tumor_vs_rest_DE_total,gene_symbol==gene)$sample_id))/length(unique(tumor_vs_rest_DE_total $sample_id))
  sample <- paste(unique(subset(tumor_vs_rest_DE_total,gene_symbol==gene)$sample_id),collapse="/")
  fc_vector <- c(fc_vector,fc)
  fdr_vector <- c(fdr_vector,fdr)
  exp_vector <- c(exp_vector,exp)
  freq_vector <- c(freq_vector,freq)
  sample_vector <- c(sample_vector,sample)
}
DE_genes_filtered_df$p_val_adj <- fdr_vector
DE_genes_filtered_df$avg_logFC <- fc_vector
DE_genes_filtered_df$avg_norm_exp <- exp_vector
DE_genes_filtered_df$sample_freq <- freq_vector
DE_genes_filtered_df$samples <- sample_vector

colnames(DE_genes_filtered_df)[1] <- "Gene"
DE_genes_filtered_surface_df <- as.data.frame(DE_genes_filtered_surface)
DE_genes_filtered_surface_df$GO_surface <- "Surface"
DE_genes_filtered_df <-  merge(DE_genes_filtered_df,DE_genes_filtered_surface_df,by.x="Gene",by.y="DE_genes_filtered_surface",all.x = TRUE)
DE_genes_filtered_df$GO_surface[is.na(DE_genes_filtered_df$GO_surface)] <- "Non-surface" 
colnames(DE_genes_filtered_df)[1] <- "Gene"

cat("SORT THE TABLE BY AVERAGE LOG FOLD CHANGE...\n")
DE_genes_filtered_df <- DE_genes_filtered_df[order(DE_genes_filtered_df$avg_logFC,decreasing = TRUE),] #sort by avg_logFC

################################################
###Cellular Loation Annotation from Cell Surface Protein atlas
################################################
#CSPA <- read.table("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/Cell_Surface_Protein_Atlas_S2_File.txt",header=TRUE,stringsAsFactors = FALSE,fill=TRUE,sep="\t")
cat("ADD CELLULAR LOCATION ANNOTATION FROM CELL SURFACE PROTEIN ATLAS AND HUMAN PROTEIN ATLAS...\n")
CSPA <- read.table(cspa_file,header=TRUE,stringsAsFactors = FALSE,fill=TRUE,sep="\t")
#CSPA_subset <- CSPA %>% filter(UniProt_Cell_surface=="yes") #CSPA_category, ENTREZ_gene_symbol
DE_genes_filtered_df$CSPA_category <- CSPA$CSPA_category[match(DE_genes_filtered_df$Gene,CSPA$ENTREZ_gene_symbol)] #very high confidence since in Uniprot cell surface

################################################
###Cellular Loation Annotation from Human Protein atlas
################################################
#HPA <- read.table("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/Human_Protein_Atlas_subcellular_location.txt",header=TRUE,stringsAsFactors = FALSE,fill=TRUE,sep="\t")
HPA <- read.table(hpa_file,header=TRUE,stringsAsFactors = FALSE,fill=TRUE,sep="\t")
HPA_subset <- filter(HPA,grepl("Plasma_membrane",Main_location))  #Gene_name,Reliability
DE_genes_filtered_df$HPA_Reliability <-  HPA_subset$Reliability[match(DE_genes_filtered_df$Gene,HPA_subset$Gene_name)]

write.table(DE_genes_filtered_df,paste0(out_path,"/",tumor_ct,"_specific_DEG_with_surface_annotations_from_3DB.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

#Gene list for GTEX tissue specificity testing
write.table(DE_genes_filtered_df$Gene,paste0(out_path,"/",tumor_ct,"_specific_DEG_with_surface_annotations_from_3DB_gene_list.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote = FALSE)

cat("RESULT IS SAVED IN ",paste0(out_path,"/",tumor_ct,"_specific_DEG_with_surface_annotations_from_3DB.txt"),"\n")
