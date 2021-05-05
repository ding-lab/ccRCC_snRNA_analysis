# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(fmsb)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")
## input the tumor clusters with low 3p loss fraction
tumorclusters_highlight_df <- fread(data.table = F, "./Resources/Analysis_Results/tumor_subcluster/assign_3ploss_status_by_tumor/20210427.v1/TumorClusters.WithLow3pLoss20210427.v1.tsv")

# preprocess for subsetting--------------------------------------------------------------
## group gene sets into modules
module1_df <- data.frame(geneset_name = c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1"),
                         module_name = "Cell_cycle")
module2_df <- data.frame(geneset_name = c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP"),
                         module_name = "Immune")
module3_df <- data.frame(geneset_name = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
                         module_name = "EMT")
module4_df <- data.frame(geneset_name = c("HALLMARK_UV_RESPONSE_DN", "HALLMARK_MTORC1_SIGNALING"),
                         module_name = "mTOR")
modules_df <- rbind(module1_df, module2_df, module3_df, module4_df)
modules_df <- modules_df %>%
  mutate(scoregroup_name =  paste0(gsub(x = geneset_name, pattern = "HALLMARK_", replacement = ""), "_Score"))
modules_df$scoregroup_name_sim <- c("Mitotic_spindle",
                                    "E2F",
                                    "G2M",
                                    "DNA_repair",
                                    "MYC",
                                    "Allograft_rejection",
                                    "Complement",
                                    "Inflammatory",
                                    "IFNγ",
                                    "KRAS",
                                    "EMT",
                                    "Hypoxia",
                                    "TNFA_signaling",
                                    "UV",
                                    "mTORC1"
                                   )
easyids_process <- unique(str_split_fixed(string = tumorclusters_highlight_df$tumor_subcluster, pattern = "_", n = 2)[,1])

# plot by sample ----------------------------------------------------------
easyid_tmp <- "C3L-00079-T1"
easyid_tmp <- "C3L-01313-T1"
for (easyid_tmp in easyids_process) {
  plotdata_pre_df <- enrich_df %>%
    mutate(cluster_name2 = gsub(x = cluster_name, pattern = "\\.", replacement = "-")) %>%
    filter(grepl(x = cluster_name2, pattern = easyid_tmp)) %>%
    mutate(cluster_label = str_split_fixed(string = cluster_name2, pattern = "_", n = 2)[,2]) %>%
    mutate(is.highlight = cluster_name2 %in% tumorclusters_highlight_df$tumor_subcluster)
  plotdata_df <- plotdata_pre_df[, modules_df$scoregroup_name]
  rownames(plotdata_df) <- plotdata_pre_df$cluster_label
  plotdata_df <- rbind(rep(250,ncol(plotdata_df)) , rep(-100,ncol(plotdata_df)) , plotdata_df)
  colnames(plotdata_df) <- modules_df$scoregroup_name_sim
  ## make colors
  colors_cluster = RColorBrewer::brewer.pal(n = nrow(plotdata_pre_df), name = "Dark2")
  names(colors_cluster) <- plotdata_pre_df$cluster_label
  ## make line types
  linetypes_cluster <- rep(1, nrow(plotdata_pre_df))
  linetypes_cluster[plotdata_pre_df$is.highlight] <- 2
  names(linetypes_cluster) <- plotdata_pre_df$cluster_label
  ## make legend label
  legend_labels <- rownames(plotdata_df[-c(1,2),])
  # legend_labels[plotdata_pre_df$is.highlight] <- paste0(legend_labels[plotdata_pre_df$is.highlight], "(without 3p loss)")
  # file2write <- paste0(dir_out, easyid_tmp, ".png")
  # png(file2write, width = 800, height = 800, res = 150)
  file2write <- paste0(dir_out, easyid_tmp, ".pdf")
  pdf(file2write, width = 5, height = 5, useDingbats = F)
  radarchart(df = plotdata_df, plty = linetypes_cluster, plwd=3,
             pcol = colors_cluster,
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
             ## title
             title = paste0(easyid_tmp, " tumor clusters")
  )
  # Add a legend
  legend(x=1, y=1.3, legend = legend_labels, 
         bty = "n", pch=20 ,  text.col = "black", cex=0.8, pt.cex = 3, col = colors_cluster)
  dev.off()
}

# plot for C3N-01200-T1 and T2 ----------------------------------------------------------
for (caseid_tmp in "C3N-01200") {
  plotdata_pre_df <- enrich_df %>%
    mutate(cluster_name2 = gsub(x = cluster_name, pattern = "\\.", replacement = "-")) %>%
    filter(grepl(x = cluster_name2, pattern = caseid_tmp)) %>%
    mutate(cluster_label = gsub(x = cluster_name2, pattern = paste0(caseid_tmp, "\\-"), replacement = "")) %>%
    mutate(is.highlight = cluster_name2 %in% tumorclusters_highlight_df$tumor_subcluster)
  plotdata_df <- plotdata_pre_df[, modules_df$scoregroup_name]
  rownames(plotdata_df) <- plotdata_pre_df$cluster_label
  plotdata_df <- rbind(rep(250,ncol(plotdata_df)) , rep(-100,ncol(plotdata_df)) , plotdata_df)
  colnames(plotdata_df) <- modules_df$scoregroup_name_sim
  ## make colors
  colors_cluster = RColorBrewer::brewer.pal(n = nrow(plotdata_pre_df), name = "Dark2")
  names(colors_cluster) <- plotdata_pre_df$cluster_label
  ## make line types
  linetypes_cluster <- rep(1, nrow(plotdata_pre_df))
  linetypes_cluster[plotdata_pre_df$is.highlight] <- 2
  names(linetypes_cluster) <- plotdata_pre_df$cluster_label
  ## make legend label
  legend_labels <- rownames(plotdata_df[-c(1,2),])
  # legend_labels[plotdata_pre_df$is.highlight] <- paste0(legend_labels[plotdata_pre_df$is.highlight], "(without 3p loss)")
  # file2write <- paste0(dir_out, caseid_tmp, ".png")
  # png(file2write, width = 800, height = 800, res = 150)
  file2write <- paste0(dir_out, caseid_tmp, ".pdf")
  pdf(file2write, width = 5, height = 5, useDingbats = F)
  radarchart(df = plotdata_df, plty = linetypes_cluster, plwd=3,
             pcol = colors_cluster, 
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
             ## title
             title = paste0(caseid_tmp, " tumor clusters")
  )
  # Add a legend
  legend(x=-1.5, y=1.3, legend = legend_labels, 
         bty = "n", pch=20 ,  text.col = "black", cex=0.7, pt.cex = 1.5, col = colors_cluster)
  dev.off()
}

