contour_hsroc <- function(x,
                          cex.axis = 1.2,
                          cex.lab = 1.2,
                          max_per_page = 3,
                          ...) {
  
  # --- Input checks ---
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hsroc")
    stop("contour_hsroc() is for HSROC models only.")
  
  samp.gen <- x$rawOutput
  if (!is.list(samp.gen))
    stop("x$rawOutput must contain a list of MCMC chains.")
  
  # --- Detect number of diagnostic tests ---
  K <- length(grep("^post\\.se\\[", colnames(samp.gen[[1]])))
  testname <- paste0("Test ", seq_len(K))
  
  if (!is.null(x$testname)) {
    testname <- x$testname
  } else if (!is.null(x$testnames)) {
    testname <- x$testnames
  } else {
    testname <- paste0("Test ", seq_len(K))
  }
  # 长度保护（防止 test 数量不一致）
  if (length(testname) != K)
    testname <- paste0("Test ", seq_len(K))
  
  plots <- vector("list", K)
  plots_per_page <- 0
  
  for (i in seq_len(K)) {
    # --- Extract posterior sensitivity/specificity ---
    se <- as.numeric(unlist(lapply(samp.gen, function(m) m[, paste0("post.se[", i, "]")])))
    sp <- as.numeric(unlist(lapply(samp.gen, function(m) m[, paste0("post.sp[", i, "]")])))
    fpr <- 1 - sp
    dat <- data.frame(x = fpr, y = se)
    
    # --- Kernel density estimation ---
    hscv1 <- tryCatch(ks::Hscv(dat), error = function(e) ks::Hpi(dat))
    fhat1 <- ks::kde(dat, H = hscv1, compute.cont = TRUE)
    dimnames(fhat1[['estimate']]) <- list(fhat1[["eval.points"]][[1]],
                                          fhat1[["eval.points"]][[2]])
    
    # --- Reshape density estimates for ggplot ---
    aa <- reshape2::melt(fhat1[['estimate']])
    aa$Var1 <- rep(fhat1$eval.points[[1]], times = length(fhat1$eval.points[[2]]))
    aa$Var2 <- rep(fhat1$eval.points[[2]], each  = length(fhat1$eval.points[[1]]))
    
    p <- ggplot(aa, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = value)) +
      geom_contour(aes(z = value), breaks = fhat1[["cont"]]["75%"]) +
      geom_contour(aes(z = value), breaks = fhat1[["cont"]]["50%"]) +
      geom_contour(aes(z = value), breaks = fhat1[["cont"]]["25%"]) +
      geom_contour(aes(z = value), breaks = fhat1[["cont"]]["10%"]) +
      geom_contour(aes(z = value), breaks = fhat1[["cont"]]["5%"]) +
      ggtitle(testname[i]) +
      xlab("False positive rate") +
      ylab("True positive rate")
    
    print(p)
    plots[[i]] <- p
    
    # --- Page handling (for multiple plots in sequence) ---
    plots_per_page <- plots_per_page + 1
    if (plots_per_page >= max_per_page && i < K) {
      grid::grid.newpage()
      plots_per_page <- 0
    }
  }
  
    invisible(plots)
}
