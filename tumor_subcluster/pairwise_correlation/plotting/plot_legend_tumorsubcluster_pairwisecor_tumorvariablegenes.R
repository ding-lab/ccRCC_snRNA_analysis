# Yige Wu @WashU Jul 2020

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


# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
### make color for methylation
methyl_color_fun <- colorRamp2(c(quantile(bulk_sn_omicsprofile_df$Methyl.VHL, 0.1, na.rm=T), 
                                 quantile(bulk_sn_omicsprofile_df$Methyl.VHL, 0.5, na.rm=T), 
                                 quantile(bulk_sn_omicsprofile_df$Methyl.VHL, 0.9, na.rm=T)),
                               c("#018571", "white", "#a6611a"))
## make colors for deletion CNV fraction
color_fun_frac_loss <-  colorRamp2(c(0, 1),c("white", cnv_state_colors["loss"]))
color_fun_frac_gain <-  colorRamp2(c(0, 1),c("white", cnv_state_colors["gain"]))
## make colors for histogical type
colors_hist_type <- c("Normal Adjacent Tissue" = "#66c2a5", "Clear cell renal cell carcinoma" = "#fc8d62", "non-Clear cell renal cell carcinoma" = "#8da0cb")
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))

# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(col_fun = col_fun, 
         title = "Pearson's coeffcient\n(variably expressed genes\nwithin tumor cells)", 
         title_position = "lefttop",
         legend_width = unit(6, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_hist_type),
         title = "Histologic Type",
         legend_gp = gpar(fill = colors_hist_type)),
  Legend(labels = c("Mutated (WES)", "Mutated (Mapped the Mutation of T1 to snRNA Reads)", "None", "No Data"),
         title = "Somatic Mutation Status",
         legend_gp = gpar(fill = c("#e7298a", "#c994c7", "white", color_na))),
  Legend(col_fun = methyl_color_fun,
         title = "Bulk VHL Promoter Methylation",
         direction = "horizontal"),
  Legend(col_fun = color_fun_frac_loss, 
         title = "Fraction of tumor cells\nwith copy number loss", 
         title_position = "lefttop",
         legend_width = unit(6, "cm"),
         direction = "horizontal"),
  Legend(col_fun = color_fun_frac_gain, 
         title = "Fraction of tumor cells\nwith copy number gain", 
         title_position = "lefttop",
         legend_width = unit(6, "cm"),
         direction = "horizontal"))

## save heatmap
file2write <- paste0(dir_out, "Legend", ".pdf")
pdf(file2write,
    width = 10, height = 5)
### combine heatmap and heatmap legend
ComplexHeatmap::draw(annotation_lgd)
dev.off()
