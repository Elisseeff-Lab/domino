#!/bin/bash
# usage example_data_03_SCENIC.sh path_to_loom_counts_data

# example_data_03_SCENIC

SCENIC_DIR="scenic"
mkdir ${SCENIC_DIR}

# install SCENIC 0.11.0 singularity image
singularity build "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" docker://aertslab/pyscenic:0.12.1

# copy input  pbmc counts data
cp $1 "${SCENIC_DIR}/"

# Matrix containing motifs as rows and genes as columns and ranking position for each gene and motif (based on CRM scores) as values. URLs provided link to v2 feather files required for the 0.12.1 version of pySENIC.
curl "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  -o "${SCENIC_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

curl "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  -o "${SCENIC_DIR}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

# List of genes encoding TFs in the HG38 reference genome
curl "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt" \
  -o "${SCENIC_DIR}/allTFs_hg38.txt"

# Motif annotations based on the 2017 cisTarget motif collection. Use these files if you are using the mc9nr databases.
curl "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
  -o "${SCENIC_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic grn \
    "${SCENIC_DIR}/pbmc3k_counts.loom" \
    "${SCENIC_DIR}/allTFs_hg38.txt" \
    -o "${SCENIC_DIR}/pbmc3k_adj.tsv" \
    --num_workers 4 \
    --seed 123

singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic ctx \
    "${SCENIC_DIR}/pbmc3k_adj.tsv" \
    "${SCENIC_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "${SCENIC_DIR}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "${SCENIC_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname "${SCENIC_DIR}/pbmc3k_counts.loom" \
    --mode "dask_multiprocessing" \
    --output "${SCENIC_DIR}/regulons_pbmc_3k.csv" \
    --num_workers 4
    
singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic aucell \
    "${SCENIC_DIR}/pbmc3k_counts.loom" \
    "${SCENIC_DIR}/regulons_pbmc_3k.csv" \
    -o "${SCENIC_DIR}/auc_pbmc_3k.csv"
