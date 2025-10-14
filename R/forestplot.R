#' @import plotrix
#' @import grDevices
#' @import plotrix
#' @import grDevices
forestplot <- function(K, nstu, samp, dat, testname, dirc) {
  index <- indicator(K, nstu, dat)
  summ  <- summary(samp)
  result <- summ[[2]]
  samp.mat <- summ[[2]]
  
  para_study_se <- paralist("stud.se", summ)
  para_study_sp <- paralist("stud.sp", summ)
  
  pool.se <- result[paste0("post.se[", 1:K, "]"), 3]
  pool.sp <- result[paste0("post.sp[", 1:K, "]"), 3]
  
  grDevices::pdf(file.path(dirc, "forest_hsroc.pdf"),
                 height = max(11, K * 5),
                 width = 7)
  oldpar <- par(mfrow = c(K, 2))
  on.exit({
    par(oldpar)
    dev.off()
  })
  
  for (k in 1:K) {
    rows_se <- para_study_se[((k - 1) * nstu + 1):(k * nstu)]
    rows_sp <- para_study_sp[((k - 1) * nstu + 1):(k * nstu)]
    
    forest.se <- samp.mat[rows_se, c(1, 3, 5), drop = FALSE]
    forest.sp <- samp.mat[rows_sp, c(1, 3, 5), drop = FALSE]
    
    # Sensitivity forest plot
    plotCI(x = forest.se[, 2], y = 1:nstu,
           li = forest.se[, 1], ui = forest.se[, 3],
           xaxt = "n", err = "x",
           xlab = paste(testname[k], " Se(%)"), ylab = "Study ID",
           slty = index[, k], xlim = c(0, 1))
    abline(v = pool.se[k], lty = 3)
    axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
    axis(2, 1:nstu, tick = TRUE)
    
    # Specificity forest plot
    plotCI(x = forest.sp[, 2], y = 1:nstu,
           li = forest.sp[, 1], ui = forest.sp[, 3],
           xaxt = "n", err = "x",
           xlab = paste(testname[k], " Sp(%)"), ylab = "Study ID",
           slty = index[, k], xlim = c(0, 1))
    abline(v = pool.sp[k], lty = 3)
    axis(1, seq(0, 1, 0.1), labels = seq(0, 100, 10), tick = TRUE)
  }
}