# Domino
Domino is a tool for analysis of intra- and intercellular signaling in single cell RNA seq based on transcription factor activation. 

##### Branch notes:
Code used in the version of Domino published in 2021 has been archived in the branch `domino_paper_archive`. Changes, issues and pull requests will be handled on the `master` branch for easy use and installation of the updated version.

##### Update notes for v0.2
Added an option to only use cells without receptor dropout when calculating Pearson correlation.

##### Citation
If you find Domino useful, please consider citing us. 

Cherry C, Maestas DR, Han J, Andorko JI, Cahan P, Fertig EJ, Garmire LX, Elisseeff JH. Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics. Nat Biomed Eng. 2021 Oct;5(10):1228-1238. doi: 10.1038/s41551-021-00770-5. Epub 2021 Aug 2. PMID: 34341534; PMCID: PMC9894531.

And make sure to cite SCENIC if you use it for your transcription factor scores.

Aibar, Sara, et al. "SCENIC: single-cell regulatory network inference and clustering." Nature methods 14.11 (2017): 1083-1086.

## Installation
Domino itself can be installed from the github repository.

    if(!require(devtools)){
        install.packages('devtools')
    }
    devtools::install_github('Chris-Cherry/domino')

Domino requires a ligand-receptor database. We use [cellphonedb](https://www.cellphonedb.org/downloads). Domino will look for four files: complexes.csv genes.csv interactions.csv and proteins.csv. If you would like to download them through terminal the commands are below.

    curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/complex_curated.csv
    mv complex_curated.csv complexes.csv
    curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/interaction_curated.csv
    mv interaction_curated.csv interactions.csv
    curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv
    mv protein_curated.csv proteins.csv
    curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv
    mv gene_input.csv genes.csv

##### Setting up SCENIC
If you choose to use SCENIC to generate transcription factor activation scores, I would recommend the Python implementation. It's significantly faster than the R version and SCENIC is the longest part of the analysis pipeline by far. Docker/Singularity images are available [here](https://pyscenic.readthedocs.io/en/latest/). 

In the tutorials, we use a bash script that we wrote to run SCENIC with their Docker image. You'll need to set up your environment properly if you wish to use our script to run SCENIC. First, you will need the docker image for SCENIC 0.10.0. The Docker image will need a directory to pull input files from and write output files to. By default, we use ~/docker_scratch, but feel free to make it whatever you like and change the directory variable in the bash script.

    mkdir ~/docker_scratch
    
Put all of the following reference files in the scratch folder if you are using our script.
    
SCENIC requires a list of [transcription factors](https://github.com/aertslab/pySCENIC/tree/master/resources), [motif annotations](https://resources.aertslab.org/cistarget/motif2tf/), and [cisTarget motifs](https://resources.aertslab.org/cistarget/) which are all available from the authors of SCENIC for human (HGNC), mouse (MGI), and fly. The following will download everything necessary for an analysis of a data set with MGI gene labels for the mm9 genome.

    curl -O https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/mm_mgi_tfs.txt
    curl -O https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-10species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-10species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-5kb-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-5kb-10species.mc9nr.feather
    
And this is for HGNC gene labels for the hg19 genome

    curl -O https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt
    curl -O https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.feather
    curl -O https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-10species.mc9nr.feather
    
## A complete basic analysis using Domino
We will run through a complete analysis with Domino picking up from phenotypic assignment of clusters. If you would like to see an example of Domino without using clusters check out the cluster free analysis vignette [here](https://github.com/Chris-Cherry/domino/vignettes). 

First we will need to get some data. [Here](https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1) is a Seurat object which contains cluster assignments from the Seurat 2,700 PBMC tutorial. You wont need to install Seurat if you don't have it already - we are just pulling data from the object. We will create a directory to put our stuff as well.

    mkdir ~/work/domino_basic_example

##### Running SCENIC
The first step is to generate transcription factor activation scores with SCENIC. To following along, you will need to have created a Docker scratch directory and downloaded the necessary HGNC reference files as described in the **Setting up SCENIC** section above. First we will need to write the counts file into the Docker temporary directory from R.

    ser = readRDS('pbmc3k_final.rds')
    write.table(t(as.matrix(ser@assays$RNA@counts)), '~/docker_scratch/counts.tsv', 
        sep = '\t', col.names = NA)
    
Now we simply need to run the SCENIC bash script for data with MGI symbols. For this data set on a 32-core EPYC Rome 7002 running this step takes 40m. If you would like to skip this step, the necessary output files are provided [here](https://www.dropbox.com/s/m4wnir4x73vpfhe/basic_analysis_scenic_outputs.zip?dl=0). Either way, move the files into your working directory by terminal or file explorer.

    bash ~/my_software/domino/scenic_bash/hgnc_scenic.sh
    mv ~/docker_scratch/auc_mtx.csv ~/work/domino_basic_example
    mv ~/docker_scratch/regulons.csv ~/work/domino_basic_example
    
##### Creating the domino object
Now that we have the SCENIC data we can run domino. For example's sake, we will pull all the necessary information from the Seurat object first but you can also simply provide the Seurat object to create_domino and the function will pull it out. We will also read in the AUC matrix from SCENIC. Note that it must be transposed from SCENIC's output.

    library(domino)
    ser = readRDS('pbmc3k_final.rds')
    counts = ser@assays$RNA@counts
    z_scores = ser@assays$RNA@scale.data
    clusters = ser@active.ident
    auc = t(read.table('auc_mtx.csv', header = TRUE, row.names = 1, 
        stringsAsFactors = FALSE, sep = ','))
    
This command generates the domino object and prepares it (compiles the LR database, calculates cluster TF enrichment, and calculates TF - R correlation). Make sure you put the correct directory for the cellphonedb data sets in.

    pbmc_dom = create_domino(signaling_db = '~/work/ref/human/cpdb/', 
        features = auc, counts = counts, z_scores = z_scores, clusters = clusters, 
        df = 'regulons.csv')
    
##### Building the signaling network and plotting options
There are four major parameters you can play with when buidling the signaling network, two for selecting TFs for each cluster and two for connecting TFs to recs. min_tf_pval is the minimum p val for a TF to be assigned to a cluster and max_tf_per_clust is the maximum number of transcription factors allowed in each cluster. The same patter is true for rec_tf_cor_threshold and max_rec_per_tf except the thresholding is on the Pearon correlation between receptor and transcription factor. Building the signaling network takes very little time so feel free to play around with these values. In general we prefer to select a pval and cor threshold such that a good portion of the TFs and recs are not being trimmed by the maximum number thresholds.

    pbmc_dom = build_domino(pbmc_dom, max_tf_per_clust = 10, 
        min_tf_pval = .001, max_rec_per_tf = 10, rec_tf_cor_threshold = .25)

Now we just need to visualize the signaling network. First, we visualize the cluster-cluster signaling with signaling_network. We are using a max threshold here (max_thresh) because the Mk cluster's expression of ITGA2B is so much higher than the other clusters/ligands it drowns other relevant signaling out. Give it a shot without the max_thresh and you'll see what I mean.

    signaling_network(pbmc_dom, edge_weight = .5, max_thresh = 2.5)
    
![Intercluster signaling network](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/intercluster_network.png)

By default, the nodes are scaled based on expression of ligands targetting them. They are larger if cells in the data set are expressing high amounts of ligands targetting the cluster. The edges are weighted based on the strength of signaling between two specific clusters. The color of the edge will match the color of the ligand cluster. For example, in the network above the teal line connecting CD8 T cells and DCs represents signaling *from* T cells *to* DCs. In order to determine the ligand-receptor-TF signaling patterns involved in that, we can zoom in on the gene network in the DCs. 

    gene_network(pbmc_dom, clust = 'DC', layout = 'fr')
    
![DC gene network](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/DC_network.png)

By default red nodes are ligands, blue are receptors, and green are transcription factors. Now lets look at the expression of the ligands that are targetting the DCs.

    incoming_signaling_heatmap(pbmc_dom, rec_clust = 'DC', max_thresh = 2.5)

![DC incoming signaling heatmap](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/DC_signaling_heatmap.png)

It looks like CCL5 is a major component of the signaling from T cells to DCs. From the gene network above it looks like the CCL5 is predicted to target LILRB1 and would activate IRF4 and IRF8 in the DCs. We can investigate some expression patterns associated with these linkages. First, lets generate a heatmap of the transcription factor activation scores.

    feat_heatmap(pbmc_dom, norm = TRUE, bool = FALSE)

![Transcription factor activation score heatmap](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/tfas_heatmap.png)

And to look at the connections between transcription factors and receptores we can create a heatmap of correlations between the two.

    cor_heatmap(pbmc_dom, bool = FALSE, mark_connections = TRUE)
    
![Rec-TF correlation heatmap](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/cor_heatmap.png)

Finally, you can visualize the entire tf-r-l network by including all clusters with gene_network. Note that the function will return the igraph object used to make the plot. This is important if you use a force directed layout and you would like to make two versions of the same plot. Here we will use default settings to allow us to see node labels and then reformat the plot with igraph options to make it less cluttered for smaller formats.

    info = gene_network(pbmc_dom, clust = levels(pbmc_dom@clusters), 
        lig_scale = FALSE, layout = 'fr')
    plot(info$graph, layout = info$layout, vertex.size = 3, edge.color = 'grey', 
        vertex.frame.color = 'black', vertex.label = NA)
    
![Default global network](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/global_network_default.png)
![Tweaked global network](https://github.com/Chris-Cherry/domino/blob/cc_working/readme_images/global_network_tweaked.png)

Thanks for taking a look at our software. If you have any questions please let us know [here](https://github.com/Chris-Cherry/domino/issues). 
    
