library(Seurat)
library(Signac)
library(lmtest)
library(future.apply)

plan("multiprocess", workers =10)
options(future.globals.maxSize = 10 * 1024^3)

###This is an RDS file with columns corresponding to peaks, and rows to barcodes
barcode=readRDS('data/Barcode2Peak.CNV.20210604.v1.RDS')

###This is an RDS file with columns corresponding to peaks, and rows to cases 
case=readRDS('data/Case2Peak.CNV.20210604.v1.RDS')

###Subset obj so we can only have Tumor and PT cells
atac=ATAC
ATAC=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% c('C3L-00004-T1', 'C3L-00010-T1',
 'C3L-00079-T1', 'C3L-00088-T2', 'C3L-00088-T1',  'C3L-00448-T1',
'C3L-00583-T1', 'C3L-00610-T1', 'C3L-00790-T1', 'C3L-00908-T1', 'C3L-00917-T1',
'C3L-01287-T1',  'C3L-01313-T1', 'C3N-00317-T1','C3N-00733-T1',
'C3N-01200-T1','C3N-01213-T1','C3L-00416-T2','C3L-01302-T1','C3L-00096-T1', 'C3L-00026-T1',
'C3N-00242-T1', 'C3N-00437-T1', 'C3N-00495-T1')) | cell_type=='PT' & Piece_ID %in%
c('C3L-00088-N','C3N-01200-N'))


####First test run is on 512 peaks only (to check, that it's working in general):

###Now check how FindMarkers modules will work:
### name it as "object"
object=ATAC
## set ident.use.1
ident.use.1='Tumor'
## set ident.use.2
ident.use.2='PT'
## make sure set Idents for the object
Idents(object)=object$cell_type

## set "features"
### first set just a few to test
features=colnames(barcode)

##assay
assay='peaksMACS2'

## input the CNV values per peak per cell
### column is named after the peak by rownames(object)
### row is named after the cells by colnames(object)
cnv_per_feature_df=barcode

# set inputs for the below process --------------------------------------------------------------
cells.1 <- WhichCells(object = object, idents = ident.use.1)
cells.2 <- WhichCells(object = object, idents = ident.use.2)
data.use=object[[assay]]
data.use <- data.use[features, c(cells.1, cells.2)]

# process inputs further --------------------------------------------------
## prepare latent.vars data frame
latent.vars <- FetchData(
  object = object,
  vars = "peak_RF_500MACS2",
  cells = c(cells.1, cells.2)
)
## prepare group.info data frame
group.info <- data.frame(row.names = c(cells.1, cells.2))
group.info[cells.1, "group"] <- "Group1"
group.info[cells.2, "group"] <- "Group2"
group.info[, "group"] <- factor(x = group.info[, "group"])

## prepare data.use object
data.use <- data.use[, rownames(group.info), drop = FALSE]
latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
latent.vars <- cbind(latent.vars, cnv_per_feature_df[rownames(group.info),])
colnames(latent.vars)=gsub('-','_',colnames(latent.vars))
rownames(data.use)=gsub('-','_',rownames(data.use))

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
    model.data <- cbind(GENE = data.use[x,], group.info, latent.vars)
    fmla <- as.formula(object = paste(
      "group ~ GENE +",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
    ))
    fmla2 <- as.formula(object = paste(
      "group ~",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
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
to.return$bonferroni=p.adjust(to.return$p_val,method='bonferroni')

to.return$chr_peak=gsub('(.*)_.*_.*','\\1',rownames(to.return))
