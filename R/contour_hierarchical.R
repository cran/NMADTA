
contour_hierarchical <- function(x,
                                 cex.axis = 1.2,
                                 cex.lab = 1.2,
                                 max_per_page = 3,
                                 ...) {
  
  # --- Input checks ---
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hierarchical")
    stop("contour_hierarchical() is for hierarchical models only.")
  
  samp.gen <- x$rawOutput
  if (!is.list(samp.gen))
    stop("x$rawOutput must contain a list of MCMC chains.")
  
  # --- Detect number of diagnostic tests ---
  cols_se <- grep("post.Se\\[", colnames(samp.gen[[1]]), value = TRUE)
  K <- length(cols_se)
  if (K == 0)
    stop("Cannot detect posterior sensitivity columns (post.Se[i]) in MCMC output.")
  
  # --- Extract test names if available ---
  if (!is.null(x$testname)) {
    testname <- x$testname
  } else if (!is.null(x$testnames)) {
    testname <- x$testnames
  } else {
    testname <- paste0("Test ", seq_len(K))
  }
  if (length(testname) != K)
    testname <- paste0("Test ", seq_len(K))
  
  plots <- vector("list", K)
  plots_per_page <- 0
  
  for (i in seq_len(K)) {
    # --- Extract posterior samples ---
    se <- as.numeric(unlist(lapply(samp.gen, function(m) m[, paste0("post.Se[", i, "]")])))
    sp <- as.numeric(unlist(lapply(samp.gen, function(m) m[, paste0("post.Sp[", i, "]")])))
    fpr <- 1 - sp
    dat <- data.frame(x = fpr, y = se)
    
    # --- Kernel density estimation ---
    hscv <- tryCatch(ks::Hscv(dat), error = function(e) ks::Hpi(dat))
    fhat <- ks::kde(dat, H = hscv, compute.cont = TRUE)
    dimnames(fhat$estimate) <- list(fhat$eval.points[[1]], fhat$eval.points[[2]])
    
    # --- Reshape density estimates ---
    aa <- reshape2::melt(fhat$estimate)
    aa$Var1 <- rep(fhat$eval.points[[1]], times = length(fhat$eval.points[[2]]))
    aa$Var2 <- rep(fhat$eval.points[[2]], each  = length(fhat$eval.points[[1]]))
    
    # --- Build contour plot (原始风格一致) ---
    p <- ggplot(aa, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = value)) +
      geom_contour(aes(z = value), breaks = fhat$cont["75%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["50%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["25%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["10%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["5%"]) +
      ggtitle(testname[i]) +
      xlab("False positive rate") +
      ylab("True positive rate")
    
    print(p)
    plots[[i]] <- p
    
    # --- Optional pagination ---
    plots_per_page <- plots_per_page + 1
    if (plots_per_page >= max_per_page && i < K) {
      grid::grid.newpage()
      plots_per_page <- 0
    }
  }
  
  invisible(plots)
}
