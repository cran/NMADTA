#' @import plotrix
#' @import grDevices
forestplot.hierarchical <- function(K, nstu, dat, testname, samp, dirc) {
  index <- indicator(K, nstu, dat)
  summ <- summary(samp)
  result <- summ[[2]]
  samp.mat <- result
  
  # Parameter names
  para_study_se <- paralist("Se.stud", summ)
  para_study_sp <- paralist("Sp.stud", summ)
  
  # Posterior pooled means
  pool.se <- result[paste0("post.Se[", seq_len(K), "]"), 3]
  pool.sp <- result[paste0("post.Sp[", seq_len(K), "]"), 3]
  
  # Open PDF
  grDevices::pdf(file.path(dirc, "forest_hierarchical.pdf"),
                 height = max(11, K * 5), width = 7)
  oldpar <- par(mfrow = c(K, 2))
  on.exit(par(oldpar))
  
  for (k in seq_len(K)) {
    idx <- ((k - 1) * nstu + 1):(k * nstu)
    forest.se <- samp.mat[para_study_se[idx], c(1, 3, 5)]
    forest.sp <- samp.mat[para_study_sp[idx], c(1, 3, 5)]
    
    # Sensitivity
    plotCI(
      x = forest.se[, 2], y = seq_len(nstu),
      li = forest.se[, 1], ui = forest.se[, 3],
      xaxt = "n", err = "x",
      xlab = paste(testname[k], " Se(%)"),
      ylab = "Study ID",
      slty = index[, k], xlim = c(0, 1)
    )
    abline(v = pool.se[k], lty = 3)
    axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
    axis(2, seq_len(nstu), tick = TRUE)
    
    # Specificity
    plotCI(
      x = forest.sp[, 2], y = seq_len(nstu),
      li = forest.sp[, 1], ui = forest.sp[, 3],
      xaxt = "n", err = "x",
      xlab = paste(testname[k], " Sp(%)"),
      ylab = "Study ID",
      slty = index[, k], xlim = c(0, 1)
    )
    abline(v = pool.sp[k], lty = 3)
    axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
  }
  
  dev.off()
}