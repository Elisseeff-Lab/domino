# Code to prepare Seurat and loom objects from pbmc3k data as example and test data
# code to prepare `pbmc` dataset
pbmc.data <- Seurat::Read10X(data.dir = "../inst/extdata/pbmc3k_filtered_gene_bc_matrices")
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
# Annotate clusters with cell phenotypes
cell_dict <- data.frame(
    cluster = c(0:8),
    cell_type = c("naive_CD4_T_cell", "CD14_monocyte", "memory_CD4_T_cell", "B_cell", "CD8_T_cell", "CD16_monocyte", "NK_cell", "dendritic_cell", "platelet")
)

pbmc$cell_type <- 
    plyr::mapvalues(
    pbmc$seurat_clusters,
    from = cell_dict$cluster,
    to = cell_dict$cell_type
    )

# save Seurat object as RDS
usethis::use_data(pbmc, overwrite = TRUE)

# save loom counts matrix
pbmc_counts <- pbmc@assays$RNA@counts
pbmc_loom <- loomR::create(filename = "../inst/extdata/pbmc3k_counts.loom", data = pbmc_counts)
pbmc_loom$close_all() # Remember to manually close connection to loom files!
