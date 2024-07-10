#' @import biomaRt
#' @import stats
#' @importFrom utils read.csv
#' @importFrom Matrix rowSums
#' @import methods
#'
NULL

#' Create a receptor - ligand map from a CellPhoneDB signaling database
#'
#' Generates a data frame of ligand-receptor interactions from a CellPhoneDB database annotating the genes encoding the interacting ligands and receptors to be queried in transcriptomic data.
#'
#' @param genes data frame or file path to table of gene names in uniprot, hgnc_symbol, or ensembl format in CellPhoneDB database format
#' @param proteins data frame or file path to table of protein features in CellPhoneDB format
#' @param interactions data frame or file path to table of protein-protein interactions in CellPhoneDB format
#' @param complexes optional: data frame or file path to table of protein complexes in CellPhoneDB format
#' @param database_name name of the database being used, stored in output
#' @param gene_conv a tuple of (from, to) or (source, target) if gene conversion to orthologs is desired; options are ENSMUSG, ENSG, MGI, or HGNC
#' @param gene_conv_host host for conversion; default ensembl, could also use mirrors if desired
#' @param alternate_convert boolean if you would like to use a non-ensembl method of conversion (must supply table; not recommended, use only if ensembl is down)
#' @param alternate_convert_table supplied table for non-ensembl method of conversion
#' @return Data frame where each row describes a possible receptor-ligand interaction
#' @export create_rl_map_cellphonedb
#' @examples
#' data(CellPhoneDB)
#' rl_map_tiny <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
#'  proteins = CellPhoneDB$proteins_tiny,
#'  interactions = CellPhoneDB$interactions_tiny,
#'  complexes =CellPhoneDB$complexes_tiny)
#' 
create_rl_map_cellphonedb <- function(
    genes, proteins, interactions, complexes = NULL, database_name = "CellPhoneDB",
    gene_conv = NULL, gene_conv_host = "https://www.ensembl.org", alternate_convert = FALSE, alternate_convert_table = NULL) {

  # Check input structures:
  check_arg(genes, c("character", "data.frame"))
  check_arg(proteins, c("character", "data.frame"))
  check_arg(interactions, c("character", "data.frame"))
  check_arg(complexes, c("character", "data.frame", "NULL"))
  check_arg(database_name, c("character"), allow_len = c(1))
  check_arg(gene_conv, c("NULL", "character"), allow_len = c(0, 2))
  check_arg(gene_conv_host, c("character"), allow_len = c(1))

  # Read in files if needed:
  genes <- read_if_char(genes)
  proteins <- read_if_char(proteins)
  interactions <- read_if_char(interactions)
  complexes <- read_if_char(complexes)

  # replace empty cells in columns annotating gene properties with 'False' There are some
  # unannotated genes in database v2.0 that seem to have been fixed in v4.0
  gene_features <- c(
    "transmembrane", "peripheral", "secreted", "secreted_highlight", "receptor",
    "integrin", "other"
  )
  proteins[proteins$receptor == "", colnames(proteins) %in% gene_features] <- "False"

  # change cases of True/False syntax from Python to TRUE/FALSE R syntax
  genes <-  conv_py_bools(genes)
  proteins <- conv_py_bools(proteins)
  interactions <- conv_py_bools(interactions)
  complexes <- conv_py_bools(complexes)

  # gene conversions
  if (!is.null(gene_conv) & !identical(gene_conv[1], gene_conv[2])) {
    # obtain conversion dictionary
    if (alternate_convert) {
      conv_dict <- table_convert_genes(genes$gene_name,
        from = gene_conv[1], to = gene_conv[2],
        alternate_convert_table
      )
    } else {
      conv_dict <- convert_genes(genes$gene_name, from = gene_conv[1], to = gene_conv[2], host = gene_conv_host)
    }
    # column 1 is the source gene names used by the reference data base column 2 is the
    # orthologous gene names for the organism to which the reference is being converted
  }
  # Step through the interactions and build rl connections.
  rl_map <- NULL
  for (i in seq_len(nrow(interactions))) {
    inter <- interactions[i, ]
    partner_a <- inter[["partner_a"]]
    partner_b <- inter[["partner_b"]]
    conversion_flag <- list()
    # features of partner_a
    a_features <- list()
    if (partner_a %in% complexes[["complex_name"]]) {
      complex_a <- complexes[complexes[["complex_name"]] == partner_a, ]
      component_a <- as.character(complex_a[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")])
      component_a <- component_a[component_a != ""]
      a_features[["uniprot_A"]] <- paste(component_a, collapse = ",")
      gene_a <- vapply(component_a, FUN.VALUE = character(1), FUN = function(x) {
        g <- unique(genes[genes[["uniprot"]] == x, c("gene_name")])
        if (!is.null(gene_conv) & !identical(gene_conv[1], gene_conv[2])) {
          # if the original gene trying to be converted is not in the gene dictionary the
          # interaction is not included in the final rl_map
          if (sum(g %in% conv_dict[, 1]) < length(g)) {
            for (gn in g) {
              conversion_flag[[gn]] <- TRUE
            }
          } else {
            g <- paste(unique(conv_dict[conv_dict[, 1] %in% g, 2]), collapse = ";")
          }
        }
        # if multiple genes are annotated for the uniprot ID, use only the first unique instance
        if(length(g) == 1){
          res <- g
        } else {
          res <- g[1]
          g_col <- paste(g, collapse = ", ")
          message(
            component_a, " has multiple encoding gene mapped in genes table.\n",
            g_col, "\n",
            "The first mapping gene is used: ", res
          )
        }
        return(res)
      })
      a_features[["gene_A"]] <- paste(gene_a, collapse = ",")
      # annotation as a receptor or ligand is based on the annotation of the complex
      a_features[["type_A"]] <- ifelse(complex_a[["receptor"]], "R", "L")
      # replace any spaces in the partner name with an underscore
      a_features[["name_A"]] <- gsub(" ", "_", partner_a)
    } else if (partner_a %in% proteins[["uniprot"]]) {
      protein_a <- proteins[proteins[["uniprot"]] == partner_a, ]
      component_a <- protein_a[["uniprot"]]
      a_features[["uniprot_A"]] <- component_a
      gene_a <- unique(genes[genes[["uniprot"]] == component_a, c("gene_name")])
      if (!is.null(gene_conv) & !identical(gene_conv[1], gene_conv[2])) {
        # if the original gene trying to be converted is not in the gene dictionary the
        # interaction is not included in the final rl_map
        if (sum(gene_a %in% conv_dict[, 1]) < length(gene_a)) {
          for (gn in gene_a) {
            conversion_flag[[gn]] <- TRUE
          }
        } else {
          gene_a <- unique(conv_dict[conv_dict[, 1] %in% gene_a, 2])
        }
      }
      gene_a <- paste(gene_a, collapse = ";")
      a_features[["gene_A"]] <- gene_a
      a_features[["type_A"]] <- ifelse(protein_a[["receptor"]], "R", "L")
      a_features[["name_A"]] <- gene_a
    } else {
      next
    }
    if (length(conversion_flag)) {
      message(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
      message(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
      next
    }
    a_df <- as.data.frame(a_features)
    # features of partner_b
    b_features <- list()
    if (partner_b %in% complexes[["complex_name"]]) {
      complex_b <- complexes[complexes[["complex_name"]] == partner_b, ]
      component_b <- as.character(complex_b[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")])
      component_b <- component_b[component_b != ""]
      b_features[["uniprot_B"]] <- paste(component_b, collapse = ",")
      gene_b <- vapply(component_b, FUN.VALUE = character(1), FUN = function(x) {
        g <- unique(genes[genes[["uniprot"]] == x, c("gene_name")])
        if (!is.null(gene_conv) & !identical(gene_conv[1], gene_conv[2])) {
          # if the original gene trying to be converted is not in the gene dictionary the
          # interaction is not included in the final rl_map
          if (sum(g %in% conv_dict[, 1]) < length(g)) {
            for (gn in g) {
              conversion_flag[[gn]] <- TRUE
            }
          } else {
            g <- paste(unique(conv_dict[conv_dict[, 1] %in% g, 2]), collapse = ";")
          }
        }
        # if multiple genes are annotated for the uniprot ID, use only the first unique instance
        if(length(g) == 1){
          res <- g
        } else {
          res <- g[1]
          g_col <- paste(g, collapse = ", ")
          message(
            component_a, " has multiple encoding gene mapped in genes table.\n",
            g_col, "\n",
            "The first mapping gene is used: ", res
          )
        }
        return(res)
      })
      b_features[["gene_B"]] <- paste(gene_b, collapse = ",")
      # annotation as a receptor or ligand is based on the annotation of the complex
      b_features[["type_B"]] <- ifelse(complex_b[["receptor"]], "R", "L")
      # replace any spaces in the partner name with an underscore
      b_features[["name_B"]] <- gsub(" ", "_", partner_b)
    } else if (partner_b %in% proteins[["uniprot"]]) {
      protein_b <- proteins[proteins[["uniprot"]] == partner_b, ]
      component_b <- protein_b[["uniprot"]]
      b_features[["uniprot_B"]] <- component_b
      gene_b <- unique(genes[genes[["uniprot"]] == component_b, c("gene_name")])
      if (!is.null(gene_conv) & !identical(gene_conv[1], gene_conv[2])) {
        # if the original gene trying to be converted is not in the gene dictionary the
        # interaction is not included in the final rl_map
        if (sum(gene_b %in% conv_dict[, 1]) < length(gene_b)) {
          for (gn in gene_b) {
            conversion_flag[[gn]] <- TRUE
          }
        } else {
          gene_b <- unique(conv_dict[conv_dict[, 1] %in% gene_b, 2])
        }
      }
      gene_a <- paste(gene_b, collapse = ";")
      b_features[["gene_B"]] <- gene_b
      b_features[["type_B"]] <- ifelse(protein_b[["receptor"]], "R", "L")
      b_features[["name_B"]] <- gene_b
    } else {
      next
    }
    if (length(conversion_flag)) {
      message(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
      message(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
      next
    }
    b_df <- as.data.frame(b_features)
    i_features <- cbind(a_df, b_df)
    i_features[["int_pair"]] <- paste(i_features[["name_A"]], i_features[["name_B"]], sep = " & ")
    i_features[["annotation_strategy"]] <- inter[["annotation_strategy"]]
    i_features[["source"]] <- inter[["source"]]
    i_features[["database_name"]] <- database_name
    rl_map <- rbind(i_features, rl_map)
  }
  # exclude rows without receptor-ligand interactions
  rl_map <- rl_map[!(rl_map$type_A == "R" & rl_map$type_B == "R") & !(rl_map$type_A == "L" & rl_map$type_B ==
    "L"), ]
  # specify column order
  rl_map <- rl_map[, c(
    "int_pair", "name_A", "uniprot_A", "gene_A", "type_A", "name_B", "uniprot_B",
    "gene_B", "type_B", "annotation_strategy", "source", "database_name"
  )]
  return(rl_map)
}

#' Create a list of genes in regulons inferred by SCENIC
#'
#' Generates a list of transcription factors and the genes targeted by the transcription factor as part of their regulon inferred by pySCENIC
#'
#' @param regulons Data frame or file path to the table of the output of the ctx function from pySCENIC
#' @return A list where names are transcription factors and the stored values are character vectors of genes in the inferred regulons
#' @export create_regulon_list_scenic
#' @examples
#' data(SCENIC)
#' regulon_list_tiny <- create_regulon_list_scenic(regulons = SCENIC$regulons_tiny)
#'
create_regulon_list_scenic <- function(regulons) {
  if (is(regulons, "character")) {
    regulons <- read.csv(regulons)
  }
  TFS <- unique(regulons[["TF"]])
  TF_targets <- lapply(TFS, function(tf) {
    regulons_small <- regulons[regulons[["TF"]] == tf, ]
    targets <- regulons_small[["TargetGenes"]]
    target_genes <- lapply(targets, function(x) {
      split_targs <- unlist(strsplit(x, ""))
      split_targs <- split_targs[seq(2, length(split_targs))]
      split_targs <- split_targs[seq(1, length(split_targs) - 1)]
      split_targs <- paste(split_targs, collapse = "")
      split_targs <- unlist(strsplit(split_targs, "['), (']"))
      split_targs_pre <- split_targs[split_targs != ""]
      split_targs_post <- split_targs_pre[seq(1, length(split_targs_pre), 2)]
      return(split_targs_post)
    })
    return(unique(unlist(target_genes)))
  })
  names(TF_targets) <- TFS
  return(TF_targets)
}

#' Create a domino object and prepare it for network construction
#'
#' This function reads in a receptor ligand signaling database, cell level
#' features of some kind (ie. output from pySCENIC), z-scored single cell data,
#' and cluster id for single cell data, calculates a correlation matrix between
#' receptors and other features (this is transcription factor module scores if
#' using pySCENIC), and finds features enriched by cluster. It will return a
#' domino object prepared for [build_domino()], which will calculate a signaling
#' network.
#'
#' @param rl_map Data frame where each row describes a receptor-ligand interaction with required columns gene_A & gene_B including the gene names for the receptor and ligand and type_A & type_B annotating if genes A and B are a ligand (L) or receptor (R)
#' @param features Either a path to a csv containing cell level features of interest (ie. the auc matrix from pySCENIC) or named matrix with cells as columns and features as rows.
#' @param counts Counts matrix for the data. This is only used to threshold receptors on dropout.
#' @param z_scores Matrix containing z-scored expression data for all cells with cells as columns and features as rows.
#' @param clusters Named factor containing cell cluster with names as cells.
#' @param use_clusters Boolean indicating whether to use clusters.
#' @param tf_targets Optional. A list where names are transcription factors and the stored values are character vectors of genes in the transcription factor's regulon.
#' @param verbose Boolean indicating whether or not to print progress during computation.
#' @param use_complexes Boolean indicating whether you wish to use receptor/ligand complexes in the receptor ligand signaling database. If FALSE, receptor/ligand pairs where either functions as a protein complex will not be considered when constructing the signaling network.
#' @param rec_min_thresh Minimum expression level of receptors by cell. Default is 0.025 or 2.5 percent of all cells in the data set. This is important when calculating correlation to connect receptors to transcription activation. If this threshold is too low then correlation calculations will proceed with very few cells with non-zero expression.
#' @param remove_rec_dropout Whether to remove receptors with 0 expression counts when calculating correlations. This can reduce false positive correlation calculations when receptors have high dropout rates.
#' @param tf_selection_method Selection of which method to target transcription factors. If 'clusters' then differential expression for clusters will be calculated. If 'variable' then the most variable transcription factors will be selected. If 'all' then all transcription factors in the feature matrix will be used. Default is 'clusters'. Note that if you wish to use clusters for intercellular signaling downstream to MUST choose clusters.
#' @param tf_variance_quantile What proportion of variable features to take if using variance to threshold features. Default is 0.5. Higher numbers will keep more features. Ignored if tf_selection_method is not 'variable'
#' @return A domino object
#' @export create_domino
#' @examples
#' example(create_rl_map_cellphonedb, echo = FALSE)
#' example(create_regulon_list_scenic, echo = FALSE)
#' data(SCENIC)
#' data(PBMC)
#'
#' pbmc_dom_tiny <- create_domino(
#'  rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
#'  counts = PBMC$RNA_count_tiny, z_scores = PBMC$RNA_zscore_tiny,
#'  clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
#'  use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE)
#'
#' pbmc_dom_tiny_no_clusters <- create_domino(
#'  rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
#'  counts = PBMC$RNA_count_tiny, z_scores =PBMC$RNA_zscore_tiny,
#'  clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
#'  use_clusters = FALSE, use_complexes = FALSE,
#'  rec_min_thresh = 0.1, remove_rec_dropout = TRUE,
#'  tf_selection_method = "all")
#'
create_domino <- function(
    rl_map, features, counts = NULL, z_scores = NULL,
    clusters = NULL, use_clusters = TRUE, tf_targets = NULL, verbose = TRUE,
    use_complexes = TRUE, rec_min_thresh = 0.025, remove_rec_dropout = TRUE,
    tf_selection_method = "clusters", tf_variance_quantile = 0.5) {

  # Check inputs:
  check_arg(rl_map, allow_class = "data.frame",
            need_vars = c("gene_A", "gene_B", "type_A", "type_B"))

  check_arg(features, allow_class = c("data.frame", "character", "matrix"))
  if (any(class(features) %in% c("data.frame", "matrix"))) {
    check_arg(features, need_rownames = TRUE, need_colnames = TRUE)
  }


  check_arg(counts, allow_class = c("matrix", "data.frame", "Matrix", "dgCMatrix"),
              need_rownames = TRUE, need_colnames = TRUE)
  check_arg(z_scores, allow_class = "matrix", need_rownames = TRUE,
              need_colnames = TRUE)
  check_arg(clusters, allow_class = "factor", need_names = TRUE)
  

  check_arg(rec_min_thresh, allow_class = c("numeric"), allow_range = c(0, 1))

  check_arg(tf_selection_method,
            allow_values = c("clusters", "variable", "all"))

  # Create object
  dom <- domino()
  dom@misc[["create"]] <- TRUE
  dom@misc[["build"]] <- FALSE
  dom@misc[["build_vars"]] <- NULL

  # Read in lr db info
  if (verbose) {
    message("Reading in and processing signaling database")
  }
  if ("database_name" %in% colnames(rl_map)) {
    dom@db_info <- rl_map
    if (verbose) {
      message("Database provided from source: ", unique(rl_map[["database_name"]]))
    }
  } else {
    dom@db_info <- rl_map
  }
  # check for receptors that match receptor complex syntax of comma seperated genes
  non_complex_index <- which(!grepl("\\,", rl_map[["gene_A"]]) & !grepl("\\,", rl_map[["gene_B"]]))
  # discard interactions including complexes if requested
  if (use_complexes == FALSE) {
    rl_map <- rl_map[non_complex_index, ]
  }
  # Get genes for receptors
  rl_reading <- NULL
  for (i in seq_len(nrow(rl_map))) {
    rl <- list()
    inter <- rl_map[i, ]
    p <- ifelse(inter[["type_A"]] == "R", "A", "B")
    q <- ifelse(p == "A", "B", "A")
    R.gene <- inter[[paste0("gene_", p)]]
    L.gene <- inter[[paste0("gene_", q)]]
    rl[["R.gene"]] <- R.gene
    rl[["L.gene"]] <- L.gene
    if (paste0("uniprot_", p) %in% names(inter)) {
      rl[["R.uniprot"]] <- inter[[paste0("uniprot_", p)]]
    }
    if (paste0("uniprot_", q) %in% names(inter)) {
      rl[["L.uniprot"]] <- inter[[paste0("uniprot_", q)]]
    }
    if (paste0("name_", p) %in% names(inter)) {
      rl[["R.name"]] <- inter[[paste0("name_", p)]]
    }
    if (paste0("name_", q) %in% names(inter)) {
      rl[["L.name"]] <- inter[[paste0("name_", q)]]
    }
    rl <- as.data.frame(rl)
    rl_reading <- rbind(rl_reading, rl)
  }
  if(nrow(rl_reading) == 0) stop("No genes annotated as receptors included in rl_map")
  # save a list of complexes and their components
  dom@linkages$complexes <- NULL
  if (use_complexes) {
    complex_list <- list()
    for (i in seq_len(nrow(rl_reading))) {
      inter <- rl_reading[i, ]
      if (grepl("\\,", inter[["L.gene"]])) {
        complex_list[[inter[["L.name"]]]] <- unlist(strsplit(inter[["L.gene"]], split = "\\,"))
      }
      if (grepl("\\,", inter[["R.gene"]])) {
        complex_list[[inter[["R.name"]]]] <- unlist(strsplit(inter[["R.gene"]], split = "\\,"))
      }
    }
    dom@linkages$complexes <- complex_list
  }
  rec_genes <- unique(unlist(strsplit(rl_reading[["R.gene"]], split = "\\,")))
  rec_names <- rl_reading[["R.name"]]
  lig_genes <- unique(unlist(strsplit(rl_reading[["L.gene"]], split = "\\,")))
  lig_names <- rl_reading[["L.name"]]
  # building RL linkages
  rec_lig_linkage <- list()
  for (rec in rec_names) {
    inter <- rl_reading[rl_reading[["R.name"]] == rec, ]
    ligs <- inter[["L.name"]]
    rec_lig_linkage[[rec]] <- ligs
  }
  dom@linkages[["rec_lig"]] <- rec_lig_linkage
  dom@misc[["rl_map"]] <- rl_reading
  # Get z-score and cluster info
  if (verbose) {
    message("Getting z_scores, clusters, and counts")
  }
  dom@z_scores <- z_scores
  if (!is.null(clusters)) {
    dom@clusters <- clusters
  }
  # Read in features matrix and calculate differential expression by cluster.
  if (is(features, "character")) {
    features <- read.csv(features, row.names = 1, check.names = FALSE)
  }
  features <- features[, colnames(dom@z_scores)]
  dom@features <- as.matrix(features)
  if (tf_selection_method == "clusters") {
    p_vals <- matrix(1, nrow = nrow(features), ncol = length(levels(dom@clusters)))
    rownames(p_vals) <- rownames(features)
    colnames(p_vals) <- levels(dom@clusters)
    if (verbose) {
      message("Calculating feature enrichment by cluster")
      clust_n <- length(levels(dom@clusters))
    }
    for (clust in levels(dom@clusters)) {
      if (verbose) {
        cur <- which(levels(dom@clusters) == clust)
        message(cur, " of ", clust_n)
      }
      cells <- which(dom@clusters == clust)
      for (feat in rownames(dom@features)) {
        p_vals[feat, clust] <- stats::wilcox.test(
          dom@features[feat, cells], dom@features[feat, -cells],
          alternative = "g"
        )$p.value
      }
    }
    dom@clust_de <- p_vals
  }
  if (tf_selection_method == "all") {
    dom@clusters <- factor()
  }
  if (tf_selection_method == "variable") {
    dom@clusters <- factor()
    variances <- apply(dom@features, 1, function(x) {
      sd(x) / mean(x)
    })
    keep_n <- length(variances) * tf_variance_quantile
    keep_id <- which(rank(variances) > keep_n)
    dom@features <- dom@features[names(keep_id), ]
  }
  # store tf_targets in linkages if they are provided as a list
  if (!is(tf_targets, "list")) {
    dom@linkages[["tf_targets"]] <- NULL
    message("tf_targets is not a list. No regulons stored")
  } else {
    dom@linkages[["tf_targets"]] <- tf_targets
  }
  # Calculate correlation matrix between features and receptors.
  dom@counts <- counts
  zero_sum <- rowSums(counts == 0)
  keeps <- which(zero_sum < (1 - rec_min_thresh) * ncol(counts))
  ser_receptors <- intersect(names(keeps), rec_genes)
  rho <- matrix(0, nrow = length(ser_receptors), ncol = nrow(dom@features))
  rownames(rho) <- ser_receptors
  colnames(rho) <- rownames(dom@features)
  if (verbose) {
    message("Calculating correlations")
    n_tf <- nrow(dom@features)
  }
  for (module in rownames(dom@features)) {
    # If df is provided then check if receptors are targets of TF. If they are then set
    # correlation equal to 0.
    if (verbose) {
      cur <- which(rownames(dom@features) == module)
      message(cur, " of ", n_tf)
    }
    if (!is.null(dom@linkages$tf_targets)) {
      tf <- gsub(pattern = "\\.\\.\\.", replacement = "", module) # correction for AUC values from pySCENIC that append an elipses to TF names due to (+) characters in the orignial python output
      module_targets <- tf_targets[[tf]]
      module_rec_targets <- intersect(module_targets, ser_receptors)
    } else {
      module_rec_targets <- NULL
    }
    scores <- dom@features[module, ]
    rhorow <- rep(0, length(ser_receptors))
    names(rhorow) <- ser_receptors
    for (rec in ser_receptors) {
      if (remove_rec_dropout) {
        keep_id <- which(dom@counts[rec, ] > 0)
        rec_z_scores <- dom@z_scores[rec, keep_id]
        tar_tf_scores <- scores[keep_id]
      } else {
        rec_z_scores <- dom@z_scores[rec, ]
        tar_tf_scores <- scores
      }
      # There are some cases where all the tfs are zero for the cells left after trimming
      # dropout for receptors. Skip those and set cor to zero manually.
      if (sum(tar_tf_scores) == 0) {
        rhorow[rec] <- 0
        next
      }
      cor <- stats::cor.test(
        rec_z_scores, tar_tf_scores, method = "spearman", 
        alternative = "greater", exact = FALSE
      )
      rhorow[rec] <- cor$estimate
    }
    if (length(module_rec_targets > 0)) {
      rhorow[module_rec_targets] <- 0
    }
    rho[, module] <- rhorow
  }
  colnames(rho) <- rownames(dom@features)
  dom@misc$rec_cor <- rho
  # assess correlation among genes in the same receptor complex
  cor_list <- list()
  for (i in seq_along(names(dom@linkages$rec_lig))) {
    r <- names(dom@linkages$rec_lig)[i]
    if (r %in% names(dom@linkages$complexes)) {
      r_genes <- dom@linkages$complexes[[r]]
    } else {
      r_genes <- r
    }
    if (sum(rownames(rho) %in% r_genes) != length(r_genes)) {
      cor_list[[r]] <- rep(0, ncol(rho))
      next
    }
    if (length(r_genes) > 1) {
      gene_cor <- rho[rownames(rho) %in% r_genes, ]
      cor_med <- apply(gene_cor, 2, function(x) {
        median(x)
      })
      cor_list[[r]] <- cor_med
    } else {
      cor_list[[r]] <- rho[rownames(rho) == r_genes, ]
    }
  }
  c_cor <- t(as.data.frame(cor_list))
  dom@cor <- c_cor
  # If cluster methods are used, calculate percentage of non-zero expression of receptor genes
  # in clusters
  if (tf_selection_method == "clusters") {
    cl_rec_percent <- NULL
    for (rec in ser_receptors) {
      rec_percent <- vapply(X = levels(dom@clusters), FUN.VALUE = numeric(1), FUN = function(x) {
        # percentage of cells in cluster with non-zero expression of receptor gene
        sum(dom@counts[rec, dom@clusters == x] > 0) / length(dom@counts[rec, dom@clusters ==
          x])
      })
      cl_rec_percent <- rbind(cl_rec_percent, rec_percent)
    }
    rownames(cl_rec_percent) <- ser_receptors
    dom@misc$cl_rec_percent <- cl_rec_percent
  }
  return(dom)
}

#' Use biomaRt to convert genes
#'
#' This function reads in a vector of genes and converts the genes to specified symbol type
#'
#' @param genes Vector of genes to convert.
#' @param from Format of gene input (ENSMUSG, ENSG, MGI, or HGNC)
#' @param to Format of gene output (MGI or HGNC)
#' @param host Host to connect to. Defaults to https://www.ensembl.org following the useMart default, but can be changed to archived hosts if useMart fails to connect.
#' @return A data frame with input genes as column 1 and converted genes as column 2
#' @keywords internal
#'
convert_genes <- function(
    genes, from = c("ENSMUSG", "ENSG", "MGI", "HGNC"), to = c("MGI", "HGNC"),
    host = "https://www.ensembl.org") {
  # Check inputs:
  stopifnot(`Genes must be a vector of characters` = (is(genes, "character") & is(genes, "vector")))
  stopifnot(`From must be one of ENSMUSG, ENSG, MGI, or HGNC` = from %in% c(
    "ENSMUSG", "ENSG", "MGI",
    "HGNC"
  ))
  stopifnot(`To must be one of MGI or HGNC` = to %in% c("MGI", "HGNC"))
  stopifnot(`Host must be  web host to connect to` = (is(host, "character") & length(host) == 1))
  if (from == "ENSMUSG") {
    srcMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    sourceAtts <- "ensembl_gene_id"
  }
  if (from == "ENSG") {
    srcMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    sourceAtts <- "ensembl_gene_id"
  }
  if (from == "MGI") {
    srcMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    sourceAtts <- "mgi_symbol"
  }
  if (from == "HGNC") {
    srcMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    sourceAtts <- "hgnc_symbol"
  }
  if (to == "MGI") {
    tarMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    tarAtts <- "mgi_symbol"
  }
  if (to == "HGNC") {
    tarMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    tarAtts <- "hgnc_symbol"
  }
  genesV2 <- getLDS(
    attributes = sourceAtts, filters = sourceAtts, values = genes, mart = srcMart,
    attributesL = tarAtts, martL = tarMart, uniqueRows = FALSE
  )
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
#' @return An updated RL signaling data frame
#' @export
#' @examples 
#' example(create_rl_map_cellphonedb, echo = FALSE)
#' lr_name <- data.frame("abbrev" = c("L", "R"), "full" = c("Ligand", "Receptor"))
#' rl_map_expanded <- add_rl_column(map = rl_map_tiny, map_ref = "type_A",
#' conv = lr_name, new_name = "type_A_full")
#' 
add_rl_column <- function(map, map_ref, conv, new_name) {
  map_in_ref <- match(map[[map_ref]], conv[, 1])
  not_in_ref <- which(is.na(map_in_ref))
  if (length(not_in_ref > 0)) {
    not_in_ref_map <- cbind.data.frame(map[not_in_ref, ], as.character(NA), stringsAsFactors = FALSE)
    colnames(not_in_ref_map)[ncol(not_in_ref_map)] <- new_name
    rownames(not_in_ref_map) <- c()
  } else {
    not_in_ref_map <- c()
  }
  new_map <- c()
  for (r_id in seq_len(nrow(map))) {
    row <- map[r_id, ]
    conv_ids <- which(conv[, 1] == row[[map_ref]])
    for (id in conv_ids) {
      new_row <- c(as.matrix(row), conv[id, 2])
      new_map <- rbind(new_map, new_row)
    }
  }
  rownames(new_map) <- c()
  colnames(new_map) <- c(colnames(map), new_name)
  new_map <- rbind.data.frame(new_map, not_in_ref_map, stringsAsFactors = FALSE)
  new_map <- data.frame(new_map, stringsAsFactors = FALSE)
}

#' Calculate mean ligand expression as a data frame for plotting in circos plot
#'
#' Creates a data frame of mean ligand expression for use in plotting a circos
#' plot of ligand expression and saving tables of mean expression.
#'us

#' @param x Gene by cell expression matrix
#' @param ligands Character vector of ligand genes to be quantified
#' @param cell_ident Vector of cell type (identity) names for which to calculate mean ligand gene expression
#' @param cell_barcodes Vector of cell barcodes (colnames of x) belonging to cell_ident to calculate mean expression across
#' @param destination Name of the receptor with which each ligand interacts
#' @return A data frame of ligand expression targeting the specified receptor
#' @export
#' @examples
#' example(build_domino, echo = FALSE)
#' counts <- dom_counts(pbmc_dom_built_tiny)
#' mean_exp <- mean_ligand_expression(counts,
#'  ligands = c("PTPRC", "FASLG"), cell_ident = "CD14_monocyte",
#'  cell_barcodes = colnames(counts), destination = "FAS")
#' 
mean_ligand_expression <- function(x, ligands, cell_ident, cell_barcodes, destination){
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
