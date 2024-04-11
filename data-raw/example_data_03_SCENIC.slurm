#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100G

module load singularity
"$(pwd)/example_data_03_SCENIC.sh" "$(pwd)/pbmc3k_counts.loom"
