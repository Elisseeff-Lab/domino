# code to prepare `pbmc3k` dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k_filtered_gene_bc_matrices")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
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
saveRDS(pbmc, file = "../data/pbmc3k_seurat.rds")

# save loom counts matrix
pbmc_counts <- pbmc@assays$RNA@counts
pbmc_loom <- loomR::create(filename = "../data/pbmc3k_counts.loom", data = pbmc_counts)
pbmc_loom$close_all() # Remember to manually close connection to loom files!
