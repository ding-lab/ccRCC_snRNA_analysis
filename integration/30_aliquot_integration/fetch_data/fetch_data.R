## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

## write table
dir_out_parent <- "/diskmnt/Projects/ccRCC_scratch/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/"
dir.create(dir_out_parent)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)
print(paste0("Output directory is at ", dir_out))

## input integrated seurat object
srat <- readRDS("/diskmnt/Projects/ccRCC_scratch/Resources/Analysis_Results/integration/30_aliquot_integration/docker_run_integration/20200211.v1/integrated_30_seurat_objects.20200211.v1.RDS")
print("Finish reading the RDS file!")
## fetch data
umap_data <- Seurat::FetchData(object = srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_data$barcode <- rownames(umap_data)
write.table(x = umap_data, file = paste0(dir_out, "30_aliquot_integration.umap_data.tsv"), quote = F, row.names = F, sep = "\t")
