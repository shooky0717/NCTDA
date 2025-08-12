Test <- function(object, X = NULL, 
                  Cell_types_to_test = NULL,
                  correction = FALSE, pv.adjust = "BY", n_neighbors = 10, order = "AMMD", 
                  n_threads = 1, cov.model = "exponential", BPPARAM = NULL,  verbose = FALSE ){
  
  if (!is.null(X)) {
    stopifnot(nrow(X) == ncol(object@gene_expression))
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }
  if(is.null(Cell_types_to_test)) {
    Cell_types_to_test <- object@cell_types
  }
  CT_to_test <- intersect(Cell_types_to_test, object@cell_types)
  CT_to_test_index <- match(intersect(Cell_types_to_test, object@cell_types), object@cell_types)
  
  if(all(is.na(CT_to_test_index))) {
    stop("No one in Cell_types_to_test matches the cell types in the object! \n")
  } else if(any(is.na(CT_to_test_index))) {
    CT_nomatch <- Cell_types_to_test[is.na(CT_to_test_index)]
    cat("The following cell types in Cell_types_to_test does not match any cell types in the object \n")
    cat(CT_nomatch, "\n")
    cat("The rest will be tested...\n")
    CT_to_test_index <- na.omit(CT_to_test_index)
  }
  
  
  y <- object@gene_expression
  X <- object@covariates
  coords <- object@location
  prop <- object@cell_type_compositions
  row_names <- rownames(y)
  n <- dim(object@gene_expression)[2] 
  K <- dim(object@cell_type_compositions)[2] 
  
  # scale coordinates proportionally
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = order, verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = n_neighbors, n_omp = 1, 
                             search.type = "tree", ordering = order_brisc, 
                             verbose = verbose)
  
  # calculate log likelihoods for nonspatial models
  nrows <- nrow(object@gene_expression)
  ncols <- ncol(object@gene_expression)
  
  loglik_lm <- vapply(seq_len(nrows), function(i) {
    y_i <- y[i, ]
    if (is.null(X)) {
      X_lm <- rep(1, ncols)
    }else{
      X_lm <- X
    }
    # model formula without intercept to enable weighted model
    as.numeric(logLik(lm(y_i ~ X_lm - 1)))
  }, numeric(1))
  
  #model fit
  
  Test1_result <- lapply(CT_to_test_index, function(i) {
    ct_test <- prop[,i]
    n <- ncols
    if(is.null(X)){
      X_0 <- matrix(1, nrow = n, ncol = 1)
    }else{
      X_0 <- X
    }
    X_full <- cbind(X_0, ct_test)
    X_full <- as.matrix(X_full)
    stopifnot(nrow(X_full) == ncol(y))
    
    ix <- seq_len(nrow(y))
    outbrisc_1 <- bplapply(ix, function(j) {
      # fit model (intercept-only model if x is NULL)
      y_j <- y[j, ]
      suppressWarnings({
        runtime <- system.time({
          out_j <- BRISC_estimation(coords = coords, y = y_j, x = X_full, 
                                    cov.model = cov.model, 
                                    ordering = order_brisc, neighbor = nn_brisc, 
                                    verbose = verbose)
        })
      })
      res_j <- c(
        out_j$Theta, 
        loglik_ut = out_j$log_likelihood, 
        runtime = runtime[["elapsed"]]
      )
      res_j
    }, BPPARAM = BPPARAM)
    matbrisc_1 <- do.call("rbind", outbrisc_1)
    # calculate statistics
    # --------------------
    
    
    matbrisc_1 <- cbind(
      matbrisc_1, 
      loglik_lm = loglik_lm
    )
    # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square
    # with 2 degrees of freedom)
    
    LR_stat <- -2 * (matbrisc_1[, "loglik_lm"] - matbrisc_1[, "loglik_ut"])
    
    pval <- 1 - pchisq(LR_stat, df = 2)
    padj <- p.adjust(pval, method = pv.adjust)
    
    
    
    output <- data.frame(matbrisc_1,
                         statistic = LR_stat,
                         pval = pval,
                         padj = padj)
    rownames(output) <- row_names
    output
  })
  names(Test1_result) <- object@cell_types[CT_to_test_index]
  ut_genes_list <- lapply(object@cell_types, function(i) {
    df <- Test1_result[[i]]
    rows_use <- df[df$padj < 0.01, ]
    rownames(rows_use)
  })
  names(ut_genes_list) <- object@cell_types
  utSVG_list <- unique(unlist(ut_genes_list))
  
  Genes_to_test <-  utSVG_list
  if(is.null(Genes_to_test)) {
    cat("All genes will be tested...\n")
    Genes_to_test <- row.names(object@gene_expression)
  }
  Genes_to_test_index <- match(intersect(Genes_to_test, row.names(object@gene_expression)), row.names(object@gene_expression))
  #if(length(Genes_to_test_index) == 1){
  #  if(is.na(Genes_to_test_index)) {
  #    stop("No one in Genes_to_test matches the genes in the object! \n")
  #  }
  #} else {
  if(all(is.na(Genes_to_test_index))) {
    stop("No one in Genes_to_test matches the genes in the object! \n")
  } else if(any(is.na(Genes_to_test_index))) {
    Genes_nomatch <- Genes_to_test[is.na(Genes_to_test_index)]
    cat("The following genes in Genes_to_test does not match any genes in the object \n")
    cat(Genes_nomatch, "\n")
    cat("The rest will be tested...\n")
    Genes_to_test_index <- na.omit(Genes_to_test_index)
  }
  
  ##fit full model
  outbrisc_2 <- bplapply(Genes_to_test_index, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[i, ]
    suppressWarnings({
      runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = X, 
                                  cov.model = cov.model, 
                                  ordering = order_brisc, neighbor = nn_brisc, 
                                  verbose = verbose)
      })
    })
    res_i <- c(
      out_i$Theta, 
      loglik_ct = out_i$log_likelihood, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = BPPARAM)
  
  # collapse output list into matrix
  matbrisc_2 <- do.call("rbind", outbrisc_2)
  
  Test2_result <- lapply(CT_to_test,function(i){
    
    loglik_ut <- Test1_result[[i]]$loglik_ut[Genes_to_test_index]
    
    # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square
    # with 2 degrees of freedom)
    
    LR_stat <- -2 * (matbrisc_2[, "loglik_ct"] - loglik_ut)
    
    pval <- 1 - pchisq(LR_stat, df = 2)
    padj <- p.adjust(pval, method = pv.adjust)
    
    # rank SVGs according to LR statistics
    LR_rank <- rank(-1 * LR_stat)
    
    output <- data.frame(matbrisc_2,
                         statistic = LR_stat,
                         rank = LR_rank,
                         pval = pval,
                         padj = padj)
    rownames(output) <-  rownames(object@gene_expression)[Genes_to_test_index]
    output
  })
  names(Test2_result) <- object@cell_types[CT_to_test_index]
  object@Test <- Test2_result
  return(object)
}
