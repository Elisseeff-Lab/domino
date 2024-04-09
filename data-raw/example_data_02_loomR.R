# example_data_02_loomR

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(loomR)

sce <- readRDS("pbmc3k_sce.rds")

# save counts data to loomR file
counts <- assay(sce, "counts")
loom <- create(
  filename = "pbmc3k_counts.loom", 
  data = as.matrix(counts)
)
loom$close_all()
