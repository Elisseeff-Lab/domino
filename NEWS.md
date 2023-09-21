# Domino2 v0.2.1

## Updates to domino object construction
- Uniform formats for inputs of receptor-ligand interaction databases, tf activity features, and regulon gene lists for operability with alternative databases and tf activity inference methods
- Helper functions for reformatting pySCENIC outputs and CellPhoneDB database files to domino-readable uniform formats
- Assessment of transcription factor linkage with receptors that function as a heteromeric complex based on correlation between transcription factor activity and all receptor component genes
- Assessment of complex ligand expression as the mean of component gene expression for plotting functions
- Minimum threshold for the percentage of cells in a cluster expressing a receptor gene for the receptor to be called active within the cluster
- Additional linkage slots for active receptors in each cluster, transcription factor-receptor linkages for each cluster, and incoming ligands for active receptors on each cluster

## Plotting functions
- Circos plot of ligand expression targeting a specified receptor where chord widths correspond to the quantity of ligand expression by each cell cluster
- Signaling networks showing only outgoing signaling from specified cell clusters
- Gene networks between two cell clusters

## Bugfixes
- Added host option for gene ortholog conversions using biomart for access to maintained mirrors
- Transcription factor-target linkages are now properly stored so that receptors in a transcription factor's regulon are excluded from linkage
- Ligand nodes sizes in gene networks correspond to quantity of ligand expression
- `create_domino()` can be run without providing a regulon list
- References to the host GitHub repository have been updated to "Elisseeff-Lab"
