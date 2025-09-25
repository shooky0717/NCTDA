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
  d <- as.matrix(dist(coords))
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = order, verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = n_neighbors, n_omp = 1, 
                             search.type = "tree", ordering = order_brisc, 
                             verbose = verbose)
  
  nrows <- nrow(object@gene_expression)
  ncols <- ncol(object@gene_expression)
  
  if (is.null(X)) {
    X_0 <- rep(1, ncols)
  }else{
    X_0 <- X
  }
  X_0tX_0 <- crossprod(X_0, X_0)
  X_0tX_0_inv <- solve(X_0tX_0 + 1e-8 * diag(ncol(X_0tX_0))) 
  P0tilde <- diag(n) - X_0 %*% X_0tX_0_inv %*% t(X_0)
  
  loglik_lm <- vapply(seq_len(nrow(y)), function(i) {
    y_i <- y[i, ]
    lm_model <- lm(y_i ~ X_0 - 1)  
    list(
      loglik = as.numeric(logLik(lm_model)),
      sigma2 = summary(lm_model)$sigma^2 
    )
  }, FUN.VALUE = list(loglik = 0, sigma2 = 0))
  loglik_lm_mat <- do.call(rbind, loglik_lm)
  rownames(loglik_lm_mat) <- row_names
  
  #model fit
  
  Test1_result <- lapply(CT_to_test_index, function(i) {
    ct_test <- prop[,i]
    n <- ncols
    X_full <- cbind(X_0, ct_test)
    X_full <- as.matrix(X_full)
    stopifnot(nrow(X_full) == ncol(y))
    
    ix <- seq_len(nrow(y))
    outbrisc_1 <- bplapply(ix, function(j) {
      y_j <- y[j, ]
      suppressWarnings({
        runtime <- system.time({
          out_j <- BRISC_estimation(coords = coords, y = y_j, x = X_full, 
                                    cov.model = cov.model, 
                                    ordering = order_brisc, neighbor = nn_brisc, 
                                    verbose = verbose)
        })
      })
      
      sigma_sq <- out_j$Theta["sigma.sq"] 
      tau_sq <- out_j$Theta["tau.sq"]     
      phi <- out_j$Theta["phi"]           
      
      if (cov.model == "exponential") {
        Sigma_spatial <- sigma_sq * exp(-d / phi)
      } else if (cov.model == "gaussian") {
        Sigma_spatial <- sigma_sq * exp(-(d / phi)^2)
      } else if (cov.model == "matern") {
        nu <- out_j$Theta["nu"]
        r <- d / phi
        Sigma_spatial <- sigma_sq * (2^(1 - nu) / gamma(nu)) * 
          (r^nu) * besselK(r, nu)
        diag(Sigma_spatial) <- sigma_sq
      }
      
      sigma2_lm_j <- loglik_lm_mat[gene_name, "sigma2"]
      P0tilde_Sigma <- P0tilde %*% Sigma_spatial
      
      e_j <- sum(diag(P0tilde_Sigma)) / sigma2_lm_j
      v_j <- sum(diag(P0tilde_Sigma %*% P0tilde_Sigma)) / (sigma2_lm_j^2)
      
      a_j <- ifelse(e_j < 1e-8, 1e-8, v_j / (2 * e_j))  
      g_j <- ifelse(v_j < 1e-8, k+1, (2 * e_j^2) / v_j)  
      
      loglik_ut_j <- out_j$log_likelihood  
      loglik_lm_j <- loglik_lm_mat[gene_name, "loglik"]  
      LR_stat_j <- -2 * (loglik_lm_j - loglik_ut_j)  
      
      pval_j <- ifelse(LR_stat_j < 0, 1, 1 - pchisq(LR_stat_j / a_j, df = g_j))
      
      c(
        loglik_ut = loglik_ut_j, loglik_lm = loglik_lm_j,
        LR_stat = LR_stat_j, sigma_sq = sigma_sq, tau_sq = tau_sq, phi = phi,
        pval = pval_j, runtime = runtime[["elapsed"]]
      )
    }, BPPARAM = BPPARAM)
    
    matbrisc_1 <- do.call(rbind, outbrisc_1)
    rownames(matbrisc_1) <- row_names
    matbrisc_1 <- cbind(Test1_df, padj = p.adjust(Test1_df[, "pval"], method = pv.adjust))
    matbrisc_1  
  })
  names(Test1_result) <- object@cell_types[CT_to_test_index]
  ut_genes_list <- lapply(object@cell_types, function(i) {
    df <- Test1_result[[i]]
    rows_use <- df[df$padj < 0.05, ]
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
  if(length(Genes_to_test_index) == 1){
   if(is.na(Genes_to_test_index)) {
    stop("No one in Genes_to_test matches the genes in the object! \n")
  }
  } else {
  if(all(is.na(Genes_to_test_index))) {
    stop("No one in Genes_to_test matches the genes in the object! \n")
  } else if(any(is.na(Genes_to_test_index))) {
    Genes_nomatch <- Genes_to_test[is.na(Genes_to_test_index)]
    cat("The following genes in Genes_to_test does not match any genes in the object \n")
    cat(Genes_nomatch, "\n")
    cat("The rest will be tested...\n")
    Genes_to_test_index <- na.omit(Genes_to_test_index)
  }
  }
  
  outbrisc_2 <- bplapply(Genes_to_test_index, function(i) {
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

  matbrisc_2 <- do.call("rbind", outbrisc_2)
  
  Test2_result <- lapply(CT_to_test,function(i){
    
    loglik_ut <- Test1_result[[i]]$loglik_ut[Genes_to_test_index]
    LR_stat <- -2 * (matbrisc_2[, "loglik_ct"] - loglik_ut)
    
    pval <- 1 - pchisq(LR_stat, df = 1)
    padj <- p.adjust(pval, method = pv.adjust)
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
