


# input dependencies ------------------------------------------------------
## input seurat object
### name it as "object"
object
## set ident.use.1
ident.use.1
## set ident.use.2
ident.use.2
## make sure set Idents for the object
Idents(object)

## set "features"
### first set just a few to test
features

## input the CNV values per peak per cell
### column is named after the peak by rownames(object)
### row is named after the cells by colnames(object)
cnv_per_feature_df

# set inputs for the below process --------------------------------------------------------------
data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
cells.1 <- WhichCells(object = object, idents = ident.use.1)
cells.2 <- WhichCells(object = object, idents = ident.use.2)
verbose = TRUE

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

# run test ----------------------------------------------------------------
my.sapply <- ifelse(
  test = verbose && nbrOfWorkers() == 1,
  yes = pbsapply,
  no = future_sapply
)
p_val <- my.sapply(
  X = rownames(x = data.use),
  FUN = function(x, cnv_df) {
    model.data <- cbind(GENE = data.use[x, ], group.info, latent.vars)
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

