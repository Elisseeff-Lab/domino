## unit tests use internal data stored in R/create_domino.R

test_that(
  "calc_rec_percentage: accurate percentage of receptor expressing cells", {
    clusters_tiny <- factor(c(rep("A", 100), rep("B", 100)))
    receptors_tiny <- c("REC1", "REC2")
    REC1 <- c(
      rep(1, 30), rep(0, 70),
      rep(1, 50), rep(0, 50)
    )
    REC2 <- c(
      rep(1, 100), rep(0, 0),
      rep(1, 10), rep(0, 90)
    )
    counts_tiny <- do.call(rbind, list(REC1, REC2))
    dimnames(counts_tiny) <- list(
      receptors_tiny,
      clusters_tiny
    )
    test_rec_percent <- calc_rec_percentage(counts = counts_tiny, clusters = clusters_tiny, receptor_genes = receptors_tiny)
    
    expect_equal(
      test_rec_percent,
      t(matrix(
        c(
          0.15, 0.25,
          0.50, 0.05
        ),
        nrow = 2, ncol = 2,
        dimnames = list(
          levels(clusters_tiny),
          receptors_tiny
        )
      ))
    )
  }
)

test_that(
  "test_tfs_rec_linkage: identify correlated TFs and receptors dataset-wide", {
    set.seed(123)
    cell_ids <- paste("cell", seq(300), sep = "_")
    clusters_tiny <- factor(
      c(rep("A", 100), rep("B", 100), rep("C", 100))
    )
    names(clusters_tiny) <- cell_ids
    simulate_poisson <- function(n, betaA, betaB, betaC) {
      id <- seq_along(3*n)
      cluster <- c(rep("A", n), rep("B", n), rep("C", n))
      E <- t(matrix(
        c(
          rep(c(1, 0, 0), n), # effect of cluster A
          rep(c(0, 1, 0), n), # effect of cluster B
          rep(c(0, 0, 1), n) # effect of cluster C
        ), 
        nrow = 3, ncol = 3 * n))
      mu <- exp(0 + betaA * E[,1] +  betaB * E[,2] +  betaC * E[,3])
      Y <- rpois(n = 3 * n, lambda = mu)
      # scale Y between values of 0 and 1
      Y_scale <- ((Y - min(Y)) / max(Y))
      df <- data.frame(
        "id" = id,
        "cluster" = cluster,
        "Y" = Y,
        "Y_scale" = Y_scale
      )
      return(df)
    }
    
    # TF1 differential increase in cluster A
    # TF2 differential increase in cluster B
    # TF3 active in all three clusters but not differential
    
    TF1 <- simulate_poisson(n = 100, betaA = 5, betaB = 1, betaC = 1)[["Y_scale"]]
    TF2 <- simulate_poisson(n = 100, betaA = 1, betaB = 5, betaC = 1)[["Y_scale"]]
    TF3 <- simulate_poisson(n = 100, betaA = 3, betaB = 3, betaC = 3)[["Y_scale"]]
    features_tiny <- t(matrix(
      c(TF1, TF2, TF3),
      ncol = 3, nrow = 300, 
      dimnames = list(
        cell_ids, 
        c("TF1", "TF2", "TF3")
      )
    ))
    feature_de_tiny <- data.frame(
      TF1 = c(0.05, 1, 1),
      TF2 = c(1, 0.05, 1),
      TF3 = c(0.5, 0.5, 0.5)
    )
    
    # REC1 is directly correlated with TF1
    # REC2 is differentially expressed in cluster B but not dependent on TF2
    # REC3 is directly correlated with TF3
    # REC4 is randomly expressed accross cells
    
    REC1 <- rpois(n = 300, lambda = TF1)
    REC2 <- simulate_poisson(n = 100, betaA = 1, betaB = 2.5, betaC = 1)[["Y"]]
    REC3 <- rpois(n = 300, lambda = TF3)
    REC4 <- rpois(n = 300, lambda = 1)
    rec_counts <- t(matrix(
      c(REC1, REC2, REC3, REC4),
      ncol = 4, nrow = 300, 
      dimnames = list(
        cell_ids, 
        c("REC1", "REC2", "REC3", "REC4")
      )
    ))
    # normalize and scale counts
    pseudo_ref <- unlist(apply((rec_counts + 1), MARGIN = 1, function(x) {exp(mean(log(x)))}))
    scale_factors <- unlist(apply((rec_counts + 1) * pseudo_ref, MARGIN = 2, median))
    rec_norm <- t(t(rec_counts) / scale_factors)
    log_rec_norm <- log(rec_norm + 1)
    # scale across genes
    rec_scale <- t(scale(t(log_rec_norm))) 
    
    rec_tf_cor <- test_tfs_rec_linkage(
      features = features_tiny, z_scores = rec_scale, counts = rec_counts, 
      feature_de = feature_de_tiny, receptors = c("REC1", "REC2", "REC3", "REC4"), 
      method = "spearman.correlation", verbose = TRUE
    )
    
    expect_true(rec_tf_cor["REC1", "TF1"] > 0.1)
    expect_true(rec_tf_cor["REC2", "TF2"] > 0.1)
    expect_true(rec_tf_cor["REC3", "TF3"] > 0.1)
    expect_true(sum(rec_tf_cor["REC4",] > 0.1) == 0)
  }
)

test_that(
  "filter_tf_regulon_receptors: set correlation statistic to 0 for receptors in regulon", { 
    cor_tiny <- diag(3) 
    dimnames(cor_tiny) <- list(
      c("REC1", "REC2", "REC3"),
      c("TF1", "TF2", "TF3")
    )
    tf_regulons_tiny <- list(
      "TF1" = c("REC1", "GENEA"),
      "TF2" = c("REC1", "GENEB"),
      "TF3" = c("REC3", "GENEC")
    )
    expected_cor <- matrix(
      c(
        0, 0, 0,
        0, 1, 0,
        0, 0, 0
      ),
      nrow = 3, ncol = 3, dimnames = list(
        c("REC1", "REC2", "REC3"),
        c("TF1", "TF2", "TF3")
      )
    )
    expect_equal(
      test_pruned_cor <- filter_tf_regulon_receptors(cor_mat = cor_tiny, tf_targets = tf_regulons_tiny), 
      expected_cor
    )
  } 
)


test_that(
  "select_cluster_tf: identify differntial feature scores by one-sided wilcox test", {
    set.seed(123)
    cell_ids <- paste("cell", seq(300), sep = "_")
    clusters_tiny <- factor(
      c(rep("A", 100), rep("B", 100), rep("C", 100))
    )
    names(clusters_tiny) <- cell_ids
    
    # there are 5 TF features
    # feature scores are between 0 and 1, generated by a poisson distribution
    # TF1 is associated with cluster A
    # TF2 is associated with cluster B
    # TF3 is associated with cluster C
    # TF4 has even elevated scores across clusters
    # TF5 has 0 score across clusters
    
    simulate_poisson <- function(n, betaA, betaB, betaC) {
      id <- seq_along(3*n)
      cluster <- c(rep("A", n), rep("B", n), rep("C", n))
      E <- t(matrix(
        c(
          rep(c(1, 0, 0), n), # effect of cluster A
          rep(c(0, 1, 0), n), # effect of cluster B
          rep(c(0, 0, 1), n) # effect of cluster C
        ), 
        nrow = 3, ncol = 3 * n))
      mu <- exp(0 + betaA * E[,1] +  betaB * E[,2] +  betaC * E[,3])
      Y <- rpois(n = 3 * n, lambda = mu)
      # scale Y between values of 0 and 1
      Y_scale <- ((Y - min(Y)) / max(Y))
      df <- data.frame(
        "id" = id,
        "cluster" = cluster,
        "Y" = Y,
        "Y_scale" = Y_scale
      )
      return(df)
    }
    
    TF1 <- simulate_poisson(n = 100, betaA = 5, betaB = 1, betaC = 1)[["Y_scale"]]
    TF2 <- simulate_poisson(n = 100, betaA = 1, betaB = 5, betaC = 1)[["Y_scale"]]
    TF3 <- simulate_poisson(n = 100, betaA = 1, betaB = 1, betaC = 5)[["Y_scale"]]
    TF4 <- simulate_poisson(n = 100, betaA = 3, betaB = 3, betaC = 3)[["Y_scale"]]
    TF5 <- simulate_poisson(n = 100, betaA = 0, betaB = 0, betaC = 0)[["Y_scale"]]
    
    features_tiny <- t(matrix(
      c(TF1, TF2, TF3, TF4, TF5),
      ncol = 5, nrow = 300, 
      dimnames = list(
        cell_ids, 
        c("TF1", "TF2", "TF3", "TF4", "TF5")
      )
    ))
    
    tf_pval_mat <- select_cluster_tf(
      features = features_tiny, clusters = clusters_tiny, 
      verbose = FALSE, method = "one.sided.wilcox"
    )
    expect_true(tf_pval_mat["TF1", "A"] < 0.05)
    expect_true(tf_pval_mat["TF2", "B"] < 0.05)
    expect_true(tf_pval_mat["TF3", "C"] < 0.05)
    expect_true(sum(tf_pval_mat["TF4",] < 0.05) == 0) 
    expect_true(sum(tf_pval_mat["TF5",] < 0.05) == 0)
  }
)

test_that(
  "write_rec_lig_linkages: generate list of receptors capable of interacting with receptors", {
    rl_parse_tiny <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC", "GENEC", "GENEB1,GENEB2", "GENEE1,GENEE2"),
      "L.gene" = c("GENEA", "GENEA", "GENED1,GENED2", "GENED1,GENED2", "GENEF"),
      "R.name" = c("complexB", "GENEC", "GENEC", "complexB", "complexE"),
      "L.name" = c("GENEA", "GENEA", "complexD", "complexD", "GENEF")
    )
    
    expected_rec_lig <- list(
      "complexB" = c("GENEA", "complexD"),
      "GENEC" = c("GENEA", "complexD"),
      "complexE" = c("GENEF")
    )
    
    expect_equal(
      write_rec_lig_linkages(rl_parse_tiny),
      expected_rec_lig
    )
  }
)
test_that(
  "read_rl_map_complexes: generate list of annotated complexes", {
    rl_parse_tiny <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC", "GENEC", "GENEB1,GENEB2"),
      "L.gene" = c("GENEA", "GENEA", "GENED1,GENED2", "GENED1,GENED2"),
      "R.name" = c("complexB", "GENEC", "GENEC", "complexB"),
      "L.name" = c("GENEA", "GENEA", "complexD", "complexD")
    )
    expected_complex_list <- list(
      "complexD" = c("GENED1", "GENED2"),
      "complexB" = c("GENEB1", "GENEB2")
    )
    expect_equal(
      read_rl_map_complexes(rl_parse_tiny, use_complexes = TRUE),
      expected_complex_list
    )
    expect_equal(
      read_rl_map_complexes(rl_parse_tiny, use_complexes = FALSE),
      NULL
    )
  }
)
test_that(
  "read_rl_map_complexes: warning that complexes with the same name have different components", {
    rl_parse_repeated <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEB1,GENEB3"),
      "L.gene" = c("GENEA", "GENEA"),
      "R.name" = c("complexB", "complexB"),
      "L.name" = c("GENEA", "GENEA")
    )
    expect_warning(read_rl_map_complexes(rl_parse_repeated, use_complexes = TRUE))
  }
)
test_that(
  "read_rl_map_genes: Convert rl_map to columns stored in dom@misc$rl_map", {
    rl_map_tiny <- data.frame(
      gene_A = c("GENEA", "GENEA"),
      protein_A = c("simpleA", "simpleA"),
      type_A = c("L", "L"),
      name_A = c("GENEA", "GENEA"),
      gene_B = c("GENEB1,GENEB2", "GENEC"),
      protein_B = c("complexB1,complexB2", "simpleC"),
      type_B = c("R", "R"),
      name_B = c("complexB", "GENEC"),
      int_pair = c("GENEA & complexB", "GENEA & GENEC"),
      annotation_strategy = c("test", "test"),
      source = c("test", "test"),
      database_name = c("CellPhoneDB", "CellPhoneDB")
    )
    expected_rl_parse <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC"),
      "L.gene" = c("GENEA", "GENEA"),
      "R.name" = c("complexB", "GENEC"),
      "L.name" = c("GENEA", "GENEA")
    )
    expect_equal(
      read_rl_map_genes(rl_map_tiny),
      expected_rl_parse
    )
  }
)
