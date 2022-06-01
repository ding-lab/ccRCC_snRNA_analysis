library(dplyr)
library(getopt)
library(optparse)

option_list = list(
  make_option(c("-o", "--output_path"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-c","--cell_type"),
              type="character",
              default=NULL,
              help="cell type of interest",
              metavar="character"),
  make_option(c("-g","--gtex"),
              type="character",
              default=NULL,
              help="tissue type of interest on GTEX",
              metavar="character"),
  make_option(c("-r","--hpaR"),
              type="character",
              default=NULL,
              help="tissue type(s) of interest on Human Protein Atlas RNA data",
              metavar="character"),
  make_option(c("-p","--hpaP"),
              type="character",
              default=NULL,
              help="tissue type of interest on Human Protein Atlas Protein data",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$output_path)){
  print_help(opt_parser)
  stop("Path to output folder is required (--output_path).n", call.=FALSE)
}
if (is.null(opt$cell_type)){
  print_help(opt_parser)
  stop("cell type of interest is required (--cell_type).n", call.=FALSE)
}
if (is.null(opt$gtex)){
  print_help(opt_parser)
  stop("tissue type of interest on GTEX data is required (--gtex).n", call.=FALSE)
}
if (is.null(opt$hpaR)){
  print_help(opt_parser)
  stop("tissue type of interest on Human Protein Atlas RNA data is required (--hpaR).n", call.=FALSE)
}
if (is.null(opt$hpaP)){
  print_help(opt_parser)
  stop("tissue type of interest on Human Protein Atlas Protein data is required (--hpaP).n", call.=FALSE)
}

# read in initial arguments
out_path <- opt$output_path
tumor_ct <- opt$cell_type
tissue_of_interest <- opt$gtex
HPA_rna_tissue_type <- strsplit(opt$hpaR, ",")[[1]]
HPA_protein_tissue_type <- opt$hpaP

dir.create(out_path,showWarnings = FALSE)
#########################################################
###GTEX tissue specificity
#########################################################
GTEX_output <- read.table(paste0(out_path,"/GTEX_tissue_specific_DEG.txt"),fill=TRUE,stringsAsFactors=F)
GTEX_filtered_genes <- GTEX_output[!grepl("---",GTEX_output$V1),]
colnames(GTEX_filtered_genes) <- GTEX_filtered_genes[1,]
GTEX_filtered_genes <- GTEX_filtered_genes[-1,c(1:5)]
colnames(GTEX_filtered_genes)[c(3,5)] <- c("PVALUE","FOLD_DIFFERENCE")
#grep captured GTEX_expression_data_subset_tumor_cells_DE_genes.txt > not_captured_GTEX_expression_data_subset_tumor_cells_DE_genes.txt
#sed -e s/#//g -i not_captured_GTEX_expression_data_subset_tumor_cells_DE_genes.txt
GTEX_cap_genes <- read.table(paste0(out_path,"/not_captured_GTEX_DEG.txt"),fill=TRUE,stringsAsFactors=F)

#be careful about below incase remove candidates accidentally
GTEX_filtered_genes <- GTEX_filtered_genes %>% filter(TISSUE==tissue_of_interest) %>% as.data.frame
write.table(GTEX_filtered_genes,paste0(out_path,"/GTEX_filtered_genes.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

comb_df <- read.table(paste0(out_path,"/",tumor_ct,"_specific_DEG_with_surface_annotations_from_3DB.txt"),sep="\t",header=T)
comb_df$GTEX_FDR <- GTEX_filtered_genes$PVALUE[match(comb_df$Gene,GTEX_filtered_genes$GENE)]
comb_df[!(comb_df$Gene %in% GTEX_cap_genes$V3) & is.na(comb_df$GTEX_FDR),"GTEX_FDR"] <- as.character("FALSE")
#write.table(comb_df,paste0(out_path,"DE_genes_filtered_surface_3DB_GTEX.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

#########################
###HPA RNA tissue specificity
########################
HPA_RNA_output <- read.table(paste0(out_path,"/HPA_RNA_tissue_specific_DEG.txt"),fill=TRUE,stringsAsFactors=F)
HPA_RNA_filtered_genes <- HPA_RNA_output[!grepl("---",HPA_RNA_output$V1),]
colnames(HPA_RNA_filtered_genes) <- HPA_RNA_filtered_genes[1,]
HPA_RNA_filtered_genes <- HPA_RNA_filtered_genes[-1,c(1:5)]
colnames(HPA_RNA_filtered_genes)[c(3,5)] <- c("PVALUE","FOLD_DIFFERENCE")

#HPA_rna_tissue_type <- c("bone_marrow","thymus","spleen","appendix","tonsil","lymph_node")
HPA_test_subset <- HPA_RNA_filtered_genes %>% filter(TISSUE %in% HPA_rna_tissue_type) %>% select(GENE,FDR)
colnames(HPA_test_subset) <- c("Gene","HPA_RNA_FDR")
HPA_test_subset <- HPA_test_subset %>% group_by(Gene) %>% summarise(HPA_RNA_FDR = min(HPA_RNA_FDR))

comb_df <- merge(comb_df,HPA_test_subset,by = "Gene",all.x = T,sort=FALSE)

#########################
###HPA tissue specificity
########################
#HPA_protein_tissue_type <- "BM_and_lymphoid_tissues"
GTEX_sig_genes <- comb_df$Gene
HPA_exp_df <- read.table("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/docker/downloaded_db/HPA_normal_tissue.tsv",header=T,fill=T)

tissue_type_match <- read.table("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/docker/downloaded_db/HPA_Tissue_type_matching.txt",sep="\t",header=T)
HPA_exp_df$general_tissue <-  tissue_type_match$general_tissue[match(HPA_exp_df$Tissue,tissue_type_match$tissue)]
HPA_exp_subset <- HPA_exp_df %>% filter(Gene_name %in% GTEX_sig_genes)

HPA_sig_genes <- read.table(paste0(out_path,"/HPA_Protein_expression_for_tissue_specific_DEG.txt"),sep="\t",fill=TRUE,header=TRUE)
HPA_sig_genes_subset <- HPA_sig_genes %>% filter(GENE_SYMBOL %in% GTEX_sig_genes)
#select the significant genes that highly expressed in lymphoid tissues
HPA_sig_genes <- subset(HPA_sig_genes,GENERAL_TISSUE_TYPE==HPA_protein_tissue_type)$GENE_SYMBOL

#find the q value for the significant genes
HPA_sig_test_df <- read.table(paste0(out_path,"/HPA_Protein_significance_test_pval_of_all_input_genes.txt"),fill=TRUE,header=TRUE)
HPA_sig_test <- subset(HPA_sig_test_df,GENE %in% GTEX_sig_genes)
HPA_non_sig_genes <- HPA_sig_test_df$GENE[!(HPA_sig_test_df$GENE %in% HPA_sig_genes)]
comb_df$HPA_Protein_FDR <- HPA_sig_test$FDR[match(comb_df$Gene,HPA_sig_test$GENE)]

#comb_df[comb_df$HPA_FDR > 0.05 & !is.na(comb_df$HPA_FDR), "HPA_FDR"] <- "FALSE"
comb_df <- comb_df %>% mutate(HPA_Protein_FDR=replace(HPA_Protein_FDR,HPA_Protein_FDR > 0.05 |  (Gene %in% HPA_non_sig_genes),"FALSE"))

write.table(comb_df,paste0(out_path,"/DE_genes_filtered_surface_3DB_GTEX_HPA.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)
cat("result is saved in", paste0(out_path,"/DE_genes_filtered_surface_3DB_GTEX_HPA.txt"),"\n")


