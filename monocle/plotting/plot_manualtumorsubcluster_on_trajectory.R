# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input paths to the monocle objects
paths_monocle_objs <- fread(data.table = F, input = "./Resources/Analysis_Results/monocle/write_paths_to_monocle_objects/20200724.v1/Paths_to_Monocle_Objects.20200724.v1.tsv")
## specify cluster number
colors_clustername <- c(RColorBrewer::brewer.pal(n = 4, name = "Dark2"), "grey50")
names(colors_clustername) <- c(paste0("C", 1:4), "NA")
swatch(colors_clustername)

# process by C3N-01200 ---------------------------------------------------------
for (id_case in "C3N-01200") {
  path_obj <- paths_monocle_objs$Path_Box[paths_monocle_objs$Case == id_case]
  ## input monocle object
  obj_monocle <- readRDS(file = path_obj)
  ## change the tumor subcluster ids to string
  pdata_df <- as.data.frame(pData(obj_monocle))
  pdata_df <- pdata_df %>%
    mutate(Name_Cluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0("C", (Id_TumorManualCluster + 1)))) %>%
    mutate(Name_AliquotCluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0(Aliquot.snRNA.WU, "_C", (Id_TumorManualCluster + 1))))
  pData(obj_monocle)$Name_Cluster <- pdata_df$Name_Cluster
  ## plot
  p <- plot_cell_trajectory(obj_monocle, color_by = "Name_Cluster",cell_size=0.6)
  p <- p + scale_color_manual(values = colors_clustername)
  p <- p + facet_wrap(~Aliquot.snRNA.WU, nrow = 2)
  p <- p + ggtitle(paste0("Pseudotime Trajectory of the Cells from Patient ", id_case))
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(aspect.ratio=1)
  ## write output
  file2write <- paste0(dir_out, id_case, ".png")
  png(file2write, width = 2000, height = 2200, res = 150)
  print(p)
  dev.off()
}

# process by C3L-00088 ---------------------------------------------------------
for (id_case in "C3L-00088") {
  path_obj <- paths_monocle_objs$Path_Box[paths_monocle_objs$Case == id_case]
  ## input monocle object
  obj_monocle <- readRDS(file = path_obj)
  ## change the tumor subcluster ids to string
  pdata_df <- as.data.frame(pData(obj_monocle))
  pdata_df <- pdata_df %>%
    mutate(Name_Cluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0("C", (Id_TumorManualCluster + 1)))) %>%
    mutate(Name_AliquotCluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0(Aliquot.snRNA.WU, "_C", (Id_TumorManualCluster + 1))))
  pData(obj_monocle)$Name_Cluster <- pdata_df$Name_Cluster
  ## plot
  p <- plot_cell_trajectory(obj_monocle, color_by = "Name_Cluster",cell_size=0.6)
  p <- p + scale_color_manual(values = colors_clustername)
  # p <- p + facet_wrap(Aliquot.snRNA.WU ~ Name_Cluster, ncol = 4)
  p <- p + facet_grid(Aliquot.snRNA.WU ~ Name_Cluster)
  p <- p + ggtitle(paste0("Pseudotime Trajectory of the Cells from Patient ", id_case))
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(aspect.ratio=1)
  ## write output
  file2write <- paste0(dir_out, id_case, ".png")
  png(file2write, width = 2000, height = 2200, res = 150)
  print(p)
  dev.off()
}

# process by C3N-01213 ---------------------------------------------------------
for (id_case in "C3N-01213") {
  path_obj <- paths_monocle_objs$Path_Box[paths_monocle_objs$Case == id_case]
  ## input monocle object
  obj_monocle <- readRDS(file = path_obj)
  ## change the tumor subcluster ids to string
  pdata_df <- as.data.frame(pData(obj_monocle))
  pdata_df <- pdata_df %>%
    mutate(Name_Cluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0("C", (Id_TumorManualCluster + 1)))) %>%
    mutate(Name_AliquotCluster = ifelse(is.na(Id_TumorManualCluster), "NA", paste0(Aliquot.snRNA.WU, "_C", (Id_TumorManualCluster + 1))))
  pData(obj_monocle)$Name_Cluster <- pdata_df$Name_Cluster
  ## plot
  p <- plot_cell_trajectory(obj_monocle, color_by = "Name_Cluster",cell_size=0.6)
  p <- p + scale_color_manual(values = colors_clustername)
  p <- p + facet_wrap(~Name_Cluster, nrow = 2)
  p <- p + ggtitle(paste0("Pseudotime Trajectory of the Cells from Patient ", id_case))
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(aspect.ratio=1)
  ## write output
  file2write <- paste0(dir_out, id_case, ".png")
  png(file2write, width = 2000, height = 2200, res = 150)
  print(p)
  dev.off()
}

