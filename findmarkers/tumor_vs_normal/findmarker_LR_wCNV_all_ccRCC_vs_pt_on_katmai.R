# Yige Wu @WashU Sep 2020
## use RNA assay according to https://github.com/satijalab/seurat/issues/2646

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(Seurat)
library(lmtest)
library(future.apply)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Data_Freezes/V2/snRNA/All_Cells_Merged/33_aliquot_merged_without_anchoring.20210428.v2.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input CNV value per barcode per gene
cnv_per_feature_df=readRDS('./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/map_CNVnex_lr_by_filteredgenes_by_snRNAbarcode/20210607.v2/Barcode2Gene.CNV.20210607.v2.RDS')

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0.1
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay <- "RNA"
DefaultAssay(srat) <- assay
cat(paste0("Assay: ", assay, "\n"))
cat("###########################################\n")
## specify test
test_process <- "LR"
## specify cell groups to compare
ident.use.1 <- "Tumor cells"
ident.use.2 <- "Proximal tubule cells from NATs"
## set "features"
### first set just a few to test
features=colnames(cnv_per_feature_df)

# preprocess the Seurat object meta data---------------------------------------------
aliquots_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Tumor" & idmetadata_df$Case != "C3L-00359"]
aliquots_tumor
aliquots_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Normal"]
aliquots_nat
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
cat("finish adding unique id for each barcode in the seurat object!\n")
## make combined id for the barcode2celltype table
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(group_findmarkers = ifelse(Cell_group5 == "Tumor cells" & orig.ident %in% aliquots_tumor, 
                                    ident.use.1, 
                                    ifelse(Cell_type.detailed == "Proximal tubule" & orig.ident %in% aliquots_nat, ident.use.2, "Others")))
## map group label
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers), warn_missing = F)
print(table(srat@meta.data$group_findmarkers))
cat("finish adding group labels\n")
Idents(srat) <- "group_findmarkers" 

# set inputs for the below process --------------------------------------------------------------
cells.1 <- WhichCells(object = srat, idents = ident.use.1)
cells.2 <- WhichCells(object = srat, idents = ident.use.2)
data.use=srat[[assay]]
data.use <- data.use[features, c(cells.1, cells.2)]

# process inputs further --------------------------------------------------
## prepare group.info data frame
group.info <- data.frame(row.names = c(cells.1, cells.2))
group.info[cells.1, "group"] <- "Group1"
group.info[cells.2, "group"] <- "Group2"
group.info[, "group"] <- factor(x = group.info[, "group"])
print(head(group.info))
cat("finish making group.info\n")

## prepare data.use object
data.use <- data.use[, rownames(group.info), drop = FALSE]
cat("finish re-ordering data.use\n")

## prepare latent.vars data frame
latent.vars <- FetchData(
  object = srat,
  vars = c("id_aliquot_barcode", "orig.ident"),
  cells = c(cells.1, cells.2)
)
latent.vars$barcode_merged <- rownames(latent.vars)
print(head(latent.vars))
# head(latent.vars$id_aliquot_barcode[!(latent.vars$id_aliquot_barcode %in% rownames(cnv_per_feature_df))])
latent.vars <- cbind(latent.vars, cnv_per_feature_df[latent.vars$id_aliquot_barcode,])
rownames(latent.vars) <- latent.vars$barcode_merged
latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
print(head(latent.vars))
cat("finish preparing latent.vars\n")

## test special symbol change
cat("printing special symbols in latent.vars\n")
colnames(latent.vars)[grepl(pattern = "\\-|\\_", x = rownames(latent.vars))]
cat("printing special symbols in data.use\n")
rownames(data.use)[grepl(pattern = "\\-|\\_", x = rownames(data.use))]
colnames(latent.vars)=gsub('-','_',colnames(latent.vars))
cat("finish changing special symbols in latent.vars\n")
rownames(data.use)=gsub('-','_',rownames(data.use))
cat("finish changing special symbols in data.use\n")

# run test ----------------------------------------------------------------
my.sapply <- ifelse(
  nbrOfWorkers() == 1,
  #  test = verbose && nbrOfWorkers() == 1,
  yes = pbsapply,
  no = future_sapply
)
p_val <- my.sapply(
  X = rownames(x = data.use),
  FUN = function(x) {
    print(x)
    model.data <- cbind(GENE = data.use[x,], group.info, latent.vars)
    fmla <- as.formula(object = paste(
      "group ~ GENE +",
      paste(c(x), collapse = "+")
    ))
    fmla2 <- as.formula(object = paste(
      "group ~",
      paste(c(x), collapse = "+")
    ))
    model1 <- glm(formula = fmla, data = model.data, family = "binomial")
    model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
    lrtest <- lrtest(model1, model2)
    return(lrtest$Pr[2])
  }
)
to.return <- data.frame(p_val, row.names = rownames(data.use))
to.return$p_val=as.numeric(as.character(unlist(to.return$p_val)))
to.return$FDR=p.adjust(to.return$p_val,method='fdr')

# write output ------------------------------------------------------------
## write output
file2write <- paste0(dir_out, test_process, 
                     ".logfc.threshold", logfc.threshold.run, 
                     ".min.pct", min.pct.run,
                     ".min.diff.pct", min.diff.pct.run,
                     ".Assay", assay,
                     ".tsv")
write.table(x = to.return, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n\n")





