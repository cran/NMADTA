forestplot.hsroc <- function(x,
                             cex.axis = 1.2, cex.lab = 1.2,
                             max_per_page = 2,
                             ...) {
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hsroc")
    stop("forestplot_hsroc() is for HSROC models only.")
  
  K        <- x$K
  nstu     <- x$nstu
  dat      <- x$dat
  testname <- x$testname
  samp     <- x$rawOutput
  
  index <- indicator(K, nstu, dat)
  summ  <- summary(samp)
  result <- summ[[2]]
  samp.mat <- summ[[2]]
  
  para_study_se <- paralist("stud.se", summ)
  para_study_sp <- paralist("stud.sp", summ)
  
  pool.se <- result[paste0("post.se[", seq_len(K), "]"), 3]
  pool.sp <- result[paste0("post.sp[", seq_len(K), "]"), 3]
  
  pages <- ceiling(K / max_per_page)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  for (p in seq_len(pages)) {
    idx_start <- (p - 1) * max_per_page + 1
    idx_end   <- min(p * max_per_page, K)
    idx <- idx_start:idx_end
    n_this <- length(idx)
    
    if (n_this > 1) {
      cex.axis <- 0.9; cex.lab <- 0.9
      margins <- c(4, 4, 2, 1) + 0.1
    } else {
      margins <- c(5, 5, 2, 2) + 0.1
    }
    
    par(mfrow = c(n_this, 2), mar = margins)
    
    for (k in idx) {
      rows_se <- para_study_se[((k - 1) * nstu + 1):(k * nstu)]
      rows_sp <- para_study_sp[((k - 1) * nstu + 1):(k * nstu)]
      
      forest.se <- samp.mat[rows_se, c(1, 3, 5), drop = FALSE]
      forest.sp <- samp.mat[rows_sp, c(1, 3, 5), drop = FALSE]
      
      plotCI(
        x = forest.se[, 2], y = seq_len(nstu),
        li = forest.se[, 1], ui = forest.se[, 3],
        xaxt = "n", err = "x",
        xlab = paste(testname[k], " Se(%)"),
        ylab = "Study ID",
        slty = index[, k], xlim = c(0, 1),
        cex.axis = cex.axis, cex.lab = cex.lab
      )
      abline(v = pool.se[k], lty = 3)
      axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
      axis(2, seq_len(nstu), tick = TRUE)
      
      plotCI(
        x = forest.sp[, 2], y = seq_len(nstu),
        li = forest.sp[, 1], ui = forest.sp[, 3],
        xaxt = "n", err = "x",
        xlab = paste(testname[k], " Sp(%)"),
        ylab = "Study ID",
        slty = index[, k], xlim = c(0, 1),
        cex.axis = cex.axis, cex.lab = cex.lab
      )
      abline(v = pool.sp[k], lty = 3)
      axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
    }
    
    if (pages > 1 && p < pages) {
      readline(prompt = sprintf("Page %d/%d done. Press [Enter] to continue...", p, pages))
    }
  }
  
  invisible(recordPlot())
}
