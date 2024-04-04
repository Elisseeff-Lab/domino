# example_data_02_loomR

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(loomR)

temp_dir <- tempdir()

# obtain sce from Zenodo
data_url <- "https://zenodo.org/records/10891532/files"

download.file(url = paste0(data_url, "/pbmc3k_sce.rds"),
              destfile = paste0(temp_dir, "/pbmc3k_sce.rds"))
sce <- readRDS(paste0(temp_dir, "/pbmc3k_sce.rds"))

# save counts data to loomR file
counts <- assay(sce, "counts")
loom <- create(
  filename = paste0(temp_dir, "/pbmc3k_counts.loom"), 
  data = as.matrix(counts)
)
loom$close_all()
