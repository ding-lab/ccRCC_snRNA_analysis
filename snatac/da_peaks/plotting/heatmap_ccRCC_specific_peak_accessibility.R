# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input accessibility data
acc_down_df=fread(data.table = F, input = './Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/Accessibility/DOWN_ccRCC_specific.Accessibility.20210617.tsv')
acc_up_df=fread(data.table = F, input ='./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/Accessibility/UP_BAP1mutants_vs_nonBAP1.Accessibility.20210612.tsv')
## input sample category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210412.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")

# make data matrix --------------------------------------------------------
## add row names
row.names(acc_down_df)=paste(acc_down_df$Group.1,acc_down_df$Group.2,sep='_')
row.names(acc_up_df)=paste(acc_up_df$Group.1,acc_up_df$Group.2,sep='_')
## make row names consistent
acc_up_df=acc_up_df[rownames(acc_down_df),]
## cbind
plotdata_df <- cbind(acc_down_df, 
                     acc_up_df[, grepl(pattern = "chr", x = colnames(acc_up_df))])
## change group names
plotdata_df$Group.1=gsub('_','-',plotdata_df$Group.1)
plotdata_df$Group.2=gsub('_','.',plotdata_df$Group.2)
## change column names
# colnames(plotdata_df)=gsub('\\.','\\-',colnames(plotdata_df))
## make ID column
plotdata_df$ID=paste(plotdata_df$Group.1,plotdata_df$Group.2,sep='_')
## change row names
rownames(plotdata_df)=plotdata_df$ID
## filter by row names
plotdata_df <- plotdata_df[!(rownames(plotdata_df) %in% c('C3L-00088-N_Tumor','C3N-01200-N_Tumor')),]
colnames(plotdata_df)[2:3]=c('Sample','Cell_type')
plotdata_df <- plotdata_df[plotdata_df$Cell_type %in% c('Tumor','PT'),]
plotdata_df2 <- plotdata_df[, colnames(plotdata_df)[grepl(pattern = "chr", x = colnames(plotdata_df))]]
plotdata_mat <- scale(x = plotdata_df2)
rownames(plotdata_mat) <- plotdata_df$Sample

# make row split ----------------------------------------------------------
row_split_vec <- plotdata_df$Cell_type
row_split_vec[row_split_vec == "Tumor"] <- "Tumor cells"
row_split_factor <- factor(plotdata_df$Cell_type, levels = c('PT','Tumor'))

# make column split -------------------------------------------------------
column_split_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_up_df), "Up peaks", "Down peaks")
column_split_factor <- factor(x = column_split_vec, levels = c("Down peaks", "Up peaks"))

# make row annotation -----------------------------------------------------
## get BAP1 mutated cases
cases_bap1 <- mut_df$Case[mut_df$mutation_category_sim %in% c("BAP1 mutated", "Both mutated")]
cases_pbrm1 <- mut_df$Case[mut_df$mutation_category_sim %in% c("PBRM1 mutated", "Both mutated")]
row_anno_df=plotdata_df %>%
  select(Sample, Cell_type) %>%
  mutate(Case = gsub(x = Sample, pattern = "\\-[A-Z][0-9]", replacement = "")) %>%
  mutate(BAP1_status = ifelse(Case %in% cases_bap1 & Cell_type == "Tumor", 'BAP1_mutant', 'NOT_BAP1_mutant')) %>%
  mutate(PBRM1_status = ifelse(Case %in% cases_pbrm1 & Cell_type == "Tumor", 'PBRM1_mutant','NOT_PBRM1_mutant'))
row_ha= rowAnnotation(#Cell_type=row_anno_df$Cell_type, 
                      BAP1_status=row_anno_df$BAP1_status,
                      PBRM1_status=row_anno_df$PBRM1_status,
                      col=list(BAP1_status=c('BAP1_mutant'='#984EA3','NOT_BAP1_mutant'='white smoke'),
                               PBRM1_status=c('PBRM1_mutant'='#FF7F00','NOT_PBRM1_mutant'='white smoke'),
                               Cell_type=c('PT'='#1B9E77','Tumor'='#E7298A')), 
                      annotation_width = unit(3, "mm"))

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.numeric(unlist(plotdata_mat)))
color_red <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[2]
colors_heatmapbody = colorRamp2(c(-2, 0, 2), 
                                c(color_blue, "white", color_red))

# plot --------------------------------------------------------------------
x=Heatmap(matrix = plotdata_mat, name='Peak\naccessibility', col = colors_heatmapbody,
          ## rows
          show_row_names = T,  show_row_dend=FALSE,  right_annotation=row_ha, row_split = row_split_factor, cluster_row_slices=F, 
          ## columns
          show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL,
          column_split = column_split_factor, cluster_column_slices=F, 
          use_raster=F)
file2write <- paste0(dir_out, "BAP1_specific_peak_accessibility.pdf")
pdf(file2write,width=6, height=5.8, useDingbats = F)
print(x)
dev.off()

x=Heatmap(matrix = plotdata_mat, name='Peak\naccessibility', col = colors_heatmapbody,
          ## rows
          show_row_names = T,  show_row_dend=FALSE,  right_annotation=row_ha, row_split = row_split_factor, cluster_row_slices=F, 
          ## columns
          show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL,
          column_split = column_split_factor, cluster_column_slices=F, 
          use_raster=T)
file2write <- paste0(dir_out, "BAP1_specific_peak_accessibility_raster.pdf")
pdf(file2write,width=6, height=5.8, useDingbats = F)
print(x)
dev.off()

