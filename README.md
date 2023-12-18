# Introducing domino2: Improved Inference of Cell Signaling from Single Cell RNA Sequencing Data <a href="https://github.com/FertigLab/domino2/tree/v0.2.1"><img src="man/figures/logo.svg" align="right" height="138" alt="domino2 repository" /></a>

domino2 is an updated version of the original [domino](https://github.com/Elisseeff-Lab/domino) R package published in Nature Biomedical Engineering in [Computational reconstruction of the signaling networks surrounding implanted biomaterials from single-cell transcriptomics](https://doi.org/10.1038/s41551-021-00770-5). domino2 is a tool for analysis of intra- and intercellular signaling in single cell RNA sequencing data based on transcription factor activation and receptor and ligand linkages.

## Installation

domino2 is undergoing active development to improve analysis capabilities and interpretability, so the codebase is subject to change as new features and fixes are implemented. **v0.2.1** of domino2 serves as the current stable version during these active updates for reproducible usage.

This version is currently hosted on the [FertigLab GitHub](https://github.com/FertigLab) as [branch v0.2.1](https://github.com/FertigLab/domino2/tree/v0.2.1) of the [domino2 repository](https://github.com/FertigLab/domino2) forked from this repository and can be installed using the remotes package.

```r
if(!require(remotes)){
    install.packages('remotes')
}
remotes::install_github('FertigLab/domino_development@v0.2.1')
```

## Usage Overview

Here is an overview of how domino2 might be used in analysis of a single cell RNA sequencing data set:

1. Transcription factor activation scores are calculated (we recommend using [pySCENIC](https://pyscenic.readthedocs.io/en/latest/), but other methods can be used as well)
2. A ligand-receptor database is used to map linkages between ligands and receptors (we recommend using [CellphoneDB](https://www.cellphonedb.org/), but other methods can be used as well).
3. A domino object is created using counts, z-scored counts, clustering information, and the data from steps 1 and 2.
4. Parameters such as the maximum number of transcription factors and receptors or the minimum correlation threshold (among others) are used to make a cell communication network
5. Communication networks can be extracted from within the domino object or visualized using a variety of plotting functions

Please see further documentation on the [domino2 repository](https://github.com/FertigLab/domino2) for additional information, such as an example analysis that includes all of these steps in detail, from downloading and running [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) to building and visualizing domino results. Other articles include further details on plotting functions and the structure of the domino object.

## Improvements
domino2 includes updates to domino object construction:
- Uniform formats for inputs of receptor-ligand interaction databases, transcription factor activity features, and regulon gene lists for operability with alternative databases and transcription factor activity inference methods
- Helper functions for reformatting pySCENIC outputs and CellPhoneDB database files to domino-readable uniform formats
- Assessment of transcription factor linkage with receptors that function as a heteromeric complex based on correlation between transcription factor activity and all receptor component genes
- Assessment of complex ligand expression as the mean of component gene expression for plotting functions
- Minimum threshold for the percentage of cells in a cluster expressing a receptor gene for the receptor to be called active within the cluster
- Additional linkage slots for active receptors in each cluster, transcription factor-receptor linkages for each cluster, and incoming ligands for active receptors on each cluster, while transcription factor-target linkages are now properly stored so that receptors in a transcription factor's regulon are excluded from linkage

Additionally, changes have been made to improve plotting functionality:
- Chord plot of ligand expression targeting a specified receptor where chord widths correspond to the quantity of ligand expression by each cell cluster
- Signaling networks showing only outgoing signaling from specified cell clusters
- Gene networks between two cell clusters
- Ligand nodes sizes in gene networks correspond to quantity of ligand expression

Some new features have been introduced:
- Addition of new class to summarize linkages in objects
- Addition of helper functions to count linkages and compare between objects
- Plotting function for differential linkages

Lastly, the package is being updated to ensure it conforms to BioConductor standards.

## Citation

If you use our package in your analysis, please cite us:

> Cherry C, Maestas DR, Han J, Andorko JI, Cahan P, Fertig EJ, Garmire LX, Elisseeff JH. Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics. Nat Biomed Eng. 2021 Oct;5(10):1228-1238. doi: 10.1038/s41551-021-00770-5. Epub 2021 Aug 2. PMID: 34341534; PMCID: PMC9894531.

> Cherry C, Mitchell J, Nagaraj S, Krishnan K, Fertig E, Elisseeff J
(2023). domino2: Cell Communication Analysis for Single Cell RNA
Sequencing. R package version 0.2.1.

## Contact Us
If you find any bugs, have questions, or want to contribute, please let us know [here](https://github.com/FertigLab/domino_development/issues).
