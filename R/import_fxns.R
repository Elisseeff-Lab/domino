#' Create a receptor-ligand map from a cellphonedb signaling database
#' 
#' DESC
#' 
#' @param genes dataframe or file path to table of gene names in uniprot, hgnc_symbol, or ensembl format in cellphonedb database format
#' @param proteins dataframe or file path to table of protein features in cellphonedb format
#' @param interactions dataframe or file path to table of protein-protein interactions in cellphonedb format
#' @param complexes optional: dataframe or file path to table of protein complexes in cellphonedb format
#' @return Data frame where each row describes a possible receptor-ligand interaction
#' @export
#' 
create_rl_map_cellphonedb = function(genes, proteins, interactions, complexes = NULL,
                                     database_name = "CellPhoneDB",
                                     gene_conv = NULL, gene_conv_host = "https://www.ensembl.org"){
  if(class(genes)[1] == "character"){
    genes = read.csv(genes, stringsAsFactors = FALSE)
  }
  if(class(proteins)[1] == "character"){
    proteins = read.csv(proteins, stringsAsFactors = FALSE)
  }
  if(class(interactions)[1] == "character"){
    interactions = read.csv(interactions, stringsAsFactors = FALSE)
  }
  if(class(complexes)[1] == "character"){
    complexes = read.csv(complexes, stringsAsFactors = FALSE)
  }
  
  # replace empty cells in columns annotating gene properties with "False"
  # There are some unannotated genes in database v2.0 that seem to have been fixed in v4.0
  gene_features = c("transmembrane", "peripheral", "secreted", "secreted_highlight", "receptor", "integrin", "other")
  proteins[proteins$receptor == "", colnames(proteins) %in% gene_features] = "False"
  
  # change cases of True/False syntax from Python to TRUE/FALSE R syntax
  for(x in colnames(genes)){
    if(identical(unique(genes[[x]]), c("True", "False")) | identical(unique(genes[[x]]), c("False", "True"))){
      genes[[x]] <- ifelse(genes[[x]] == "True", TRUE, FALSE)
    }
  }
  for(x in colnames(proteins)){
    if(identical(unique(proteins[[x]]), c("True", "False")) | identical(unique(proteins[[x]]), c("False", "True"))){
      proteins[[x]] <- ifelse(proteins[[x]] == "True", TRUE, FALSE)
    }
  }
  for(x in colnames(interactions)){
    if(identical(unique(interactions[[x]]), c("True", "False")) | identical(unique(interactions[[x]]), c("False", "True"))){
      interactions[[x]] <- ifelse(interactions[[x]] == "True", TRUE, FALSE)
    }
  }
  if(!is.null(complexes)){
    for(x in colnames(complexes)){
      if(identical(unique(complexes[[x]]), c("True", "False")) | identical(unique(complexes[[x]]), c("False", "True"))){
        complexes[[x]] <- ifelse(complexes[[x]] == "True", TRUE, FALSE)
      }
    }
  }
  
  # gene conversions
  if(!is.null(gene_conv)){
    # obtain conversion dictionary
    conv_dict = convert_genes(
      genes$gene_name, from = gene_conv[1], to = gene_conv[2], host = gene_conv_host)
    # column 1 is the source gene names used by the reference data base
    # column 2 is the orthologous gene names for the organism to which the reference is being converted
  }
  
  # Step through the interactions and build rl connections.
  rl_map <- NULL
  
  for(i in 1:nrow(interactions)){
    inter = interactions[i,]
    partner_a = inter[["partner_a"]]
    partner_b = inter[["partner_b"]]
    conversion_flag = list()
    
    # features of partner_a
    a_features <- list()
    if(partner_a %in% complexes[["complex_name"]]){
      complex_a = complexes[complexes[["complex_name"]] == partner_a,]
      component_a = as.character(
        complex_a[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
      )
      component_a = component_a[component_a != ""]
      a_features[["uniprot_A"]] = paste(component_a, collapse = ",")
      gene_a = sapply(component_a, function(x){
        g = unique(genes[genes[["uniprot"]] == x, c("gene_name")])
        if(!is.null(gene_conv)){
          # if the original gene trying to be converted is not in the gene dictionary
          # the interaction is not included in the final rl_map
          if(sum(g %in% conv_dict[,1]) < length(g)){
            for(gn in g) {conversion_flag[[gn]] = TRUE}
          }
          else{g = paste(conv_dict[conv_dict[,1] %in% g, 2], collapse = ";")}
        }
        return(g)
        }
      )
      a_features[["gene_A"]] = paste(gene_a, collapse = ",")
      # annotation as a receptor or ligand is based on the annotation of the complex
      a_features[["type_A"]] = ifelse(complex_a[["receptor"]], "R", "L")
      a_features[["name_A"]] = partner_a
    } else if(partner_a %in% proteins[["uniprot"]]) {
      protein_a = proteins[proteins[["uniprot"]] == partner_a,]
      component_a = protein_a[["uniprot"]]
      a_features[["uniprot_A"]] = component_a
      gene_a = unique(genes[genes[["uniprot"]] == component_a, c("gene_name")])
      if(!is.null(gene_conv)){
        # if the original gene trying to be converted is not in the gene dictionary
        # the interaction is not included in the final rl_map
        if(sum(gene_a %in% conv_dict[,1]) < length(gene_a)){
          for(gn in gene_a) {conversion_flag[[gn]] = TRUE}
        }
        else{gene_a = conv_dict[conv_dict[,1] %in% gene_a, 2]}
      }
      gene_a = paste(gene_a, collapse = ";")
      a_features[["gene_A"]] = gene_a
      a_features[["type_A"]] = ifelse(protein_a[["receptor"]], "R", "L")
      a_features[["name_A"]] = gene_a
    } else {
      next
    }
    if(length(conversion_flag)){
      print(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
      print(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
      next
    }
    a_df <- as.data.frame(a_features)
    
    # features of partner_b
    b_features <- list()
    if(partner_b %in% complexes[["complex_name"]]){
      complex_b = complexes[complexes[["complex_name"]] == partner_b,]
      component_b = as.character(
        complex_b[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
      )
      component_b = component_b[component_b != ""]
      b_features[["uniprot_B"]] = paste(component_b, collapse = ",")
      gene_b = sapply(component_b, function(x){
        g = unique(genes[genes[["uniprot"]] == x, c("gene_name")])
        if(!is.null(gene_conv)){
          # if the original gene trying to be converted is not in the gene dictionary
          # the interaction is not included in the final rl_map
          if(sum(g %in% conv_dict[,1]) < length(g)){
            for(gn in g) {conversion_flag[[gn]] = TRUE}
          }
          else{g = paste(conv_dict[conv_dict[,1] %in% g, 2], collapse = ";")}
        }
        return(g)
      }
      )
      b_features[["gene_B"]] = paste(gene_b, collapse = ",")
      # annotation as a receptor or ligand is based on the annotation of the complex
      b_features[["type_B"]] = ifelse(complex_b[["receptor"]], "R", "L")
      b_features[["name_B"]] = partner_b
    } else if(partner_b %in% proteins[["uniprot"]]) {
      protein_b = proteins[proteins[["uniprot"]] == partner_b,]
      component_b = protein_b[["uniprot"]]
      b_features[["uniprot_B"]] = component_b
      gene_b = unique(genes[genes[["uniprot"]] == component_b, c("gene_name")])
      if(!is.null(gene_conv)){
        # if the original gene trying to be converted is not in the gene dictionary
        # the interaction is not included in the final rl_map
        if(sum(gene_b %in% conv_dict[,1]) < length(gene_b)){
          for(gn in gene_b) {conversion_flag[[gn]] = TRUE}
        }
        else{gene_b = conv_dict[conv_dict[,1] %in% gene_b, 2]}
      }
      gene_a = paste(gene_b, collapse = ";")
      b_features[["gene_B"]] = gene_b
      b_features[["type_B"]] = ifelse(protein_b[["receptor"]], "R", "L")
      b_features[["name_B"]] = gene_b
    } else {
      next
    }
    if(length(conversion_flag)){
      print(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
      print(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
      next
    }
    b_df = as.data.frame(b_features)
    i_features = cbind(a_df, b_df)
    
    i_features[["int_pair"]] = 
      paste(i_features[["name_A"]], i_features[["name_B"]], sep = " & ")
    i_features[["annotation_strategy"]] = inter[["annotation_strategy"]]
    i_features[["source"]] = inter[["source"]]
    i_features[["database_name"]] = database_name
    rl_map <- rbind(i_features, rl_map)
  }
  # exclude rows without receptor-ligand interactions
  rl_map <- rl_map[!(rl_map$type_A == "R" & rl_map$type_B == "R") &
                     !(rl_map$type_A == "L" & rl_map$type_B == "L"),]
  
  # specify column order
  rl_map <- rl_map[, c("int_pair", 
                       "name_A", "uniprot_A", "gene_A", "type_A",
                       "name_B", "uniprot_B", "gene_B", "type_B",
                       "annotation_strategy", "source", "database_name")]
  return(rl_map)
}


#' Create a domino object and prepare it for network construction
#' 
#' This function reads in a receptor ligand signaling database, cell level 
#' features of some kind (ie. output from pySCENIC), z-scored single cell data, 
#' and cluster id for single cell data, calculates a correlation matrix between 
#' receptors and other features (this is transcription factor module scores if 
#' using pySCENIC), and finds features enriched by cluster. It will return a 
#' domino object prepared for build_domino, which will calculate a signaling 
#' network.
#' 
#' @param signaling_db Path to directory of signaling database directory. The directory must include genes.csv, proteins.csv, interactions.csv, and complexes.csv formated according to cellphonedb2 syntax.
#' @param features Either a path to a csv containing cell level features of interest (ie. the auc matrix from pySCENIC) or named matrix with cells as columns and features as rows.
#' @param ser A Seurat object containing scaled RNA expression data in the RNA assay slot and cluster identity. Either a ser object OR z_scores and clusters must be provided. If ser is present z_scores and clusters will be ignored.
#' @param counts The counts matrix for the data. If a Seurat object is provided this will be ignored. This is only used to threshold receptors on dropout.
#' @param z_scores A matrix containing z-scored expression data for all cells with cells as columns and features as rows. Either z_scores and clusters must be provided OR a ser object. If ser is present z_scores and clusters will be ignored.
#' @param clusters A named factor containing cell cluster with names as cells. Either clusters and z_scores OR ser must be provided. If ser is present z_scores and clusters will be ignored.
#' @param use_clusters Boolean indicating whether to use the clusters from a Seurat object. If a Seurat object is not provided then this parameter is ignored.
#' @param df Optional. Either a path to discovered motifs from pySCENIC as a csv file or a data frame following the format of df.csv from pySCENIC
#' @param gene_conv Optional. Vector of length two containing some combination of 'ENSMUSG', 'ENSG', 'MGI', or 'HGNC' where the first vector is the current gene format in the database and the second is the gene format in the data set. If present, the function will use biomaRt to convert the database to the data sets gene format.
#' @param gene_conv_host Optional. Host to connect to when using gene_conv. Defaults to https://www.ensembl.org following the useMart default, but can be changed to archived hosts if useMart fails to connect.
#' @param verbose Boolean indicating whether or not to print progress during computation.
#' @param use_complexes Boolean indicating whether you wish to use receptor/ligand complexes in the receptor ligand signaling database. This may lead to problems if genes which are preserved acrossed many functionally different signaling complexes are found highly expressed or correlated with features in your data set.
#' @param rec_min_thresh Minimum expression level of receptors by cell. Default is 0.025 or 2.5 percent of all cells in the data set. This is important when calculating correlation to connect receptors to transcription activation. If this threshold is too low then correlation calculations will proceed with very few cells with non-zero expression.
#' @param remove_rec_dropout Whether to remove receptors with 0 expression counts when calculating correlations. This can reduce false positive correlation calculations when receptors have high dropout rates.
#' @param tf_selection_method Selection of which method to target transcription factors. If 'clusters' then differential expression for clusters will be calculated. If 'variable' then the most variable transcription factors will be selected. If 'all' then all transcription factors in the feature matrix will be used. Default is 'clusters'. Note that if you wish to use clusters for intercellular signaling downstream to MUST choose clusters.
#' @param tf_variance_quantile What proportion of variable features to take if using variance to threshold features. Default is 0.5. Higher numbers will keep more features. Ignored if tf_selection_method is not 'variable'
#' @return A domino object.
#' @export
#'
create_domino = function(signaling_db, features, ser = NULL, counts = NULL, 
    z_scores = NULL, clusters = NULL, use_clusters = TRUE, df = NULL, 
    gene_conv = NULL, gene_conv_host = "https://www.ensembl.org",
    verbose = TRUE, use_complexes = TRUE, 
    rec_min_thresh = .025, remove_rec_dropout = TRUE, 
    tf_selection_method = 'clusters', tf_variance_quantile = .5){

    dom = domino()
    dom@misc[['tar_lr_cols']] = c('R.orig', 'L.orig')
    dom@misc[['create']] = TRUE
    dom@misc[['build']] = FALSE
    dom@misc[['build_vars']] = NULL
    if(!is.null(ser) & (!is.null(clusters) | !is.null(z_scores) | !is.null(counts))){
        warning("Ser and z_score, clusters, or counts provided. Defaulting to ser.")
    }
    if(is.null(ser) & (is.null(clusters) | is.null(z_scores) | is.null(counts))){
        stop("Either ser or clusters and z_scores must be provided")
    }
    if(!(tf_selection_method %in% c('all', 'clusters', 'variable'))){
        stop("tf_selection_method must be one of all, clusters, or variable")
    }

    # Read in lr db info
    if(verbose){print('Reading in and processing signaling database')}
    genes = read.csv(paste0(signaling_db, '/genes.csv'),
        stringsAsFactors = FALSE)
    proteins = read.csv(paste0(signaling_db, '/proteins.csv'),
        stringsAsFactors = FALSE)
    interactions = read.csv(paste0(signaling_db, '/interactions.csv'),
        stringsAsFactors = FALSE)
    complexes = read.csv(paste0(signaling_db, '/complexes.csv'),
        stringsAsFactors = FALSE)
    
    dom@db_info = list(genes = genes, proteins = proteins, 
        interactions = interactions, complexes = complexes)

    # Get all possible rec/ligs
    rec_uniprot = proteins$uniprot[which(proteins$receptor == 'True')]
    lig_uniprot = proteins$uniprot[which(proteins$receptor != 'True')]
    if(use_complexes){
        rec_complexes = complexes$complex_name[
            which(complexes$receptor == TRUE)]
        lig_complexes = complexes$complex_name[
            which(complexes$receptor == FALSE)]
    }
    
    # Step through the interactions and build rl connections.
    rl_map = matrix(0, ncol = 2, nrow = 0)
    colnames(rl_map) = c('R.uniprot', 'L.uniprot')
    for(i in 1:nrow(interactions)){
        a_pname = interactions$protein_name_a[i]
        b_pname = interactions$protein_name_b[i]

        # If not using complexes, make sure both members are proteins or next
        if(!use_complexes){
            if(a_pname == '' | b_pname == ''){
                next
            }
        }
        aid = interactions$partner_a[i]
        bid = interactions$partner_b[i]
        # Figure out recs and ligs
        # Partner A
        # See if it's receptor/ligand complex/protein and add to rec/ligs
        if(length(which(rec_uniprot == aid)) != 0){
            recs = rec_uniprot[which(rec_uniprot == aid)]
        } else if(length(which(lig_uniprot == aid)) != 0){
            ligs = lig_uniprot[which(lig_uniprot == aid)]
        } else if(length(which(rec_complexes == aid)) != 0){
            comp = rec_complexes[which(rec_complexes == aid)]
            recs = complexes[which(complexes[,1] == comp), c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
            recs = recs[-which(recs == '')]
        } else if(length(which(lig_complexes == aid)) != 0){
            comp = lig_complexes[which(lig_complexes == aid)]
            ligs = complexes[which(complexes[,1] == comp), c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
            ligs = ligs[-which(ligs == '')]
        } else {
            stop(paste('Partner A has no comp or prot match in row', i))
        }

        # Partner B
        if(length(which(rec_uniprot == bid)) != 0){
            recs = rec_uniprot[which(rec_uniprot == bid)]
        } else if(length(which(lig_uniprot == bid)) != 0){
            ligs = lig_uniprot[which(lig_uniprot == bid)]
        } else if(length(which(rec_complexes == bid)) != 0){
            comp = rec_complexes[which(rec_complexes == bid)]
            recs = complexes[which(complexes[,1] == comp), c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
            recs = recs[-which(recs == '')]
        } else if(length(which(lig_complexes == bid)) != 0){
            comp = lig_complexes[which(lig_complexes == bid)]
            ligs = complexes[which(complexes[,1] == comp), c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
            ligs = ligs[-which(ligs == '')]
        } else {
            stop(paste('Partner B has no comp or prot match in row', i))
        }

        # Add them all to the map
        for(l in ligs){
            for(r in recs){
                rl_map = rbind(rl_map, c(r, l))
            }
        }
    }


    # Get genes for receptors
    rec_prots = as.character(unique(rl_map[,1]))
    rec_orig = genes[match(rec_prots, genes[,2]),3]
    rl_map = data.frame(rl_map)
    rl_map = add_rl_column(rl_map, 'R.uniprot', cbind(rec_prots, rec_orig), 
        'R.orig')

    # Get genes for ligands
    lig_prots = as.character(unique(rl_map[,2]))
    lig_orig = genes[match(lig_prots, genes[,2]),3]
    rl_map = add_rl_column(rl_map, 'L.uniprot', cbind(lig_prots, lig_orig), 
        'L.orig')

    # Convert if needed
    if(!is.null(gene_conv)){
        conv = convert_genes(rl_map$R.orig, from = gene_conv[1], to = gene_conv[2], 
                             host = gene_conv_host)
        rl_map = add_rl_column(rl_map, 'R.orig', conv, 'R.conv')
        conv = convert_genes(rl_map$L.orig, from = gene_conv[1], to = gene_conv[2], 
                             host = gene_conv_host)
        rl_map = add_rl_column(rl_map, 'L.orig', conv, 'L.conv')
        dom@misc[['tar_lr_cols']] = c('R.conv', 'L.conv')
    }

    # Remove duplicate rows then build rl linkage
    rl_map = unique.data.frame(rl_map)
    rl_list = list()
    tar_map = rl_map[, dom@misc$tar_lr_cols]
    for(rec in unique(tar_map[,1])){
        if(is.na(rec)){next}
        ligs = unique(tar_map[which(tar_map[,1] == rec),2])
        for(lig in ligs){
            if(!is.na(lig)){
                rl_list[[rec]] = c(rl_list[[rec]], lig)
            }
        }
    }
    dom@linkages[['rec_lig']] = rl_list
    dom@misc[['rl_map']] = rl_map

    # Get z-score and cluster info
    if(verbose){print('Getting z_scores, clusters, and counts')}
    if(!is.null(ser)){
        z_scores = ser@assays$RNA@scale.data
        if(use_clusters){
            clusters = ser@active.ident
        }
        counts = ser@assays$RNA@counts
    }

    dom@z_scores = z_scores
    if(!is.null(clusters)){
        dom@clusters = clusters
    }

    # Read in features matrix and calculate differential expression by cluster.
    if(class(features)[1] == 'character'){
        features = read.csv(features, row.names = 1, check.names = FALSE)
    }
    features = features[, colnames(dom@z_scores)]
    dom@features = as.matrix(features)

    if(tf_selection_method == 'clusters'){
        p_vals = matrix(1.0, nrow = nrow(features), 
            ncol = length(levels(dom@clusters)))
        rownames(p_vals) = rownames(features)
        colnames(p_vals) = levels(dom@clusters)

        if(verbose){
            print('Calculating feature enrichment by cluster')
            clust_n = length(levels(dom@clusters))
        }
        for(clust in levels(dom@clusters)){
            if(verbose){
                cur = which(levels(dom@clusters) == clust)
                print(paste0(cur, ' of ', clust_n))
            }
            cells = which(dom@clusters == clust)
            for(feat in rownames(dom@features)){
                p_vals[feat, clust] = wilcox.test(dom@features[feat, cells], 
                    dom@features[feat, -cells], alternative = 'g')$p.value
            }
        }

        dom@clust_de = p_vals
    }

    if(tf_selection_method == 'all'){
        dom@clusters = factor()
    }

    if(tf_selection_method == 'variable'){
        dom@clusters = factor()
        variances = apply(dom@features, 1, function(x){
            sd(x)/mean(x)
        })
        keep_n = length(variances) * tf_variance_quantile
        keep_id = which(rank(variances) > keep_n)
        dom@features = dom@features[names(keep_id),]
    }

    # If present, read in and process df
    if(class(df)[1] == 'character'){
        df = read.csv(df, skip = 3, header = FALSE, stringsAsFactors = FALSE)
    }
    if(!is.null(df)){
        tf_targets = list()
        for(row in 1:nrow(df)){
            tf = df[row,1]
            tf_genes = df[row,10]
            tf_genes = gsub("[\\(']", "", 
                regmatches(tf_genes, gregexpr("\\('.*?'", tf_genes))[[1]])
            tf_targets[[tf]] = unique(c(tf_targets[[tf]], tf_genes))
        }
        dom@linkages[['tf_targets']] = tf_targets
    } else {
        dom@linkages[['tf_targets']] = NULL
    }

    # Calculate correlation matrix between features and receptors.
    dom@counts = counts
    all_receptors = unique(names(dom@linkages$rec_lig))
    zero_sum = Matrix::rowSums(counts == 0)
    keeps = which(zero_sum < .975*ncol(counts))
    ser_receptors = intersect(names(keeps), all_receptors)
    rho = matrix(0, nrow = length(ser_receptors), ncol = nrow(dom@features))
    rownames(rho) = ser_receptors
    colnames(rho) = rownames(dom@features)
    if(verbose){
        print('Calculating correlations')
        n_tf = nrow(dom@features)
    }
    for(module in rownames(dom@features)){
        # If df is provided then check if receptors are targets of TF. If they
        # are then set correlation equal to 0.
        if(verbose){
            cur = which(rownames(dom@features) == module)
            print(paste0(cur, ' of ', n_tf))
        }
        if(!is.null(dom@linkages$tf_targets)){
            tf = substring(module, first = 1, last = nchar(module) - 3)
            module_targets = tf_targets[[tf]]
            module_rec_targets = intersect(module_targets, ser_receptors)
        }
        scores = dom@features[module,]
        rhorow = rep(0, length(ser_receptors))
        names(rhorow) = ser_receptors
        for(rec in ser_receptors){
            if(remove_rec_dropout){
                keep_id = which(dom@counts[rec,] > 0)
                rec_z_scores = dom@z_scores[rec,keep_id]
                tar_tf_scores = scores[keep_id]
            } else {
                rec_z_scores = dom@z_scores[rec,]
                tar_tf_scores = scores
            }

            # There are some cases where all the tfas are zero for the cells
            # left after trimming dropout for receptors. Skip those and set
            # cor to zero manually.
            if(sum(tar_tf_scores) == 0){
                rhorow[rec] = 0
                next
            }
            cor = cor.test(rec_z_scores, tar_tf_scores, 
                method = 'spearman', alternative = 'greater')
            rhorow[rec] = cor$estimate
        }
        if(length(module_rec_targets > 0)){
            rhorow[module_rec_targets] = 0
        }
        rho[,module] = rhorow
    }
    colnames(rho) = rownames(dom@features)
    dom@cor = rho
    return(dom)
}

#' Use biomaRt to convert genes
#' 
#' This function reads in a vector of genes and converts the genes to 
#' 
#' @param genes Vector of genes to convert.
#' @param from Format of gene input (ENSMUSG, ENSG, MGI, or HGNC)
#' @param to Format of gene output (ENSMUSG, ENSG, MGI, or HGNC)
#' @param host Host to connect to. Defaults to https://www.ensembl.org following the useMart default, but can be changed to archived hosts if useMart fails to connect.
#' @return A data frame with input genes as col 1 and output as col 2.
#' 
convert_genes = function(genes, from, to, host = "https://www.ensembl.org"){
    if (from == 'ENSMUSG'){
        srcMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                                   host = host)
        sourceAtts =  'ensembl_gene_id'   
    }
    if (from == 'ENSG'){
        srcMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                                   host = host)
        sourceAtts = 'ensembl_gene_id'
    }
    if (from == 'MGI'){
        srcMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                                   host = host)
        sourceAtts = 'mgi_symbol'    
    }
    if (from == 'HGNC'){
        srcMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                                   host = host)
        sourceAtts = 'hgnc_symbol'
    }
    if (to == 'MGI'){
        tarMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                                   host = host)
        tarAtts = 'mgi_symbol'
    }
    if (to == 'HGNC'){
        tarMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                                   host = host)
        tarAtts = 'hgnc_symbol'
    }
    genesV2 = biomaRt::getLDS(attributes = sourceAtts, filters = sourceAtts,
                     values = genes, mart = srcMart, 
                     attributesL = tarAtts, martL = tarMart,
                     uniqueRows = F)
    return(genesV2)
}

#' Adds a column to the RL signaling data frame.
#' 
#' This function adds a column to the internal rl 'map' used to map all 
#' receptor and receptor complexes to all ligand and ligand complexes.
#' 
#' @param map RL signaling data frame.
#' @param map_ref Name of column to match new data to
#' @param conv Data frame matching current data in map to new data.
#' @param new_name Name of new column to be created in RL map
#' @return RL signaling data frame.
#'
add_rl_column = function(map, map_ref, conv, new_name){
    map_in_ref = match(map[[map_ref]], conv[,1])
    not_in_ref = which(is.na(map_in_ref))
    if(length(not_in_ref > 0)){
        not_in_ref_map = cbind.data.frame(map[not_in_ref,], as.character(NA), 
            stringsAsFactors = FALSE)
        colnames(not_in_ref_map)[ncol(not_in_ref_map)] = new_name
        rownames(not_in_ref_map) = c()
    } else {
        not_in_ref_map = c()
    }
    new_map = c()
    for(r_id in 1:nrow(map)){
        row = map[r_id,]
        conv_ids = which(conv[,1] == row[[map_ref]])
        for(id in conv_ids){
            new_row = c(as.matrix(row), conv[id, 2])
            new_map = rbind(new_map, new_row)
        }
    }
    rownames(new_map) = c()
    colnames(new_map) = c(colnames(map), new_name)
    new_map = rbind.data.frame(new_map, not_in_ref_map, stringsAsFactors = FALSE)
    new_map = data.frame(new_map, stringsAsFactors = FALSE)
}

#' Calculate mean ligand expression as a data.frame for plotting in circos plot
#' 
#' Creates a data frame of mean ligand expression for use in plotting a circos
#' plot of ligand exression and saving tables of mean expression.
#' 
#' @param x gene by cell expression matrix
#' @param ligands character vector of ligand genes to be quantified
#' @param cell_ident
#' @param cell_barcodes vector of cell barcodes (colnames of x) belonging to cell_ident to calculate mean expression across
#' @param destination name of the receptor with which each ligand interacts
#' @return data frame of ligand expression targeting the specified receptor
#'
mean_ligand_expression <- 
  function(x, ligands, cell_ident, cell_barcodes, destination){
    # initiate data frame to store results
    df <- NULL
    
    for(feat in ligands){
      # index of ligand row
      lig_index <- grep(paste0("^", feat, "$"), rownames(x))
      # column indecies of cells belonging to cell_ident
      cell_index <- colnames(x) %in% cell_barcodes
      
      cell_df <- data.frame(
        origin = paste0(cell_ident, "_", feat),
        destination = destination,
        mean.expression = mean(x[lig_index, cell_index])
      )
      
      df <- rbind(df, cell_df)
    }
    return(df)
  }

