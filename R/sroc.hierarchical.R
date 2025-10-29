
sroc.hierarchical <- function(x,
                              cex.axis = 1.2,
                              cex.lab  = 1.2,
                              max_per_page = 3,
                              ...) {
  # --- Input checks ---
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hierarchical")
    stop("sroc.hierarchical() is for hierarchical models only.")
  
  # --- Extract pieces ---
  samp <- x$rawOutput
  dat  <- x$dat
  nstu <- x$nstu
  K    <- x$K
  if (!is.list(samp))
    stop("x$rawOutput must contain a list of MCMC chains.")
  
  # --- test names ---
  if (!is.null(x$testname)) {
    testname <- x$testname
  } else if (!is.null(x$testnames)) {
    testname <- x$testnames
  } else {
    testname <- paste0("Test ", seq_len(K))
  }
  if (length(testname) != K) testname <- paste0("Test ", seq_len(K))
  
  # --- Core setup ---
  index <- indicator(K, nstu, dat)
  l      <- except(index, K, nstu)
  summ   <- summary(samp)
  result <- summ[[2]]
  para_study_se <- paralist("Se.stud", summ)
  para_study_sp <- paralist("Sp.stud", summ)
  stud.se <- matrix(result[para_study_se, 3], ncol = K, nrow = nstu, byrow = FALSE)
  stud.sp <- matrix(result[para_study_sp, 3], ncol = K, nrow = nstu, byrow = FALSE)
  smry <- summ$statistics[, c("Mean", "SD")]
  cov.id <- grep("Cov", rownames(smry))
  cov <- smry[cov.id, "Mean"]
  
  # --- helpers ---
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  safe_keep <- function(lj, n) {
    bad <- unlist(lj, use.names = FALSE)
    if (length(bad) == 0) return(seq_len(n))
    setdiff(seq_len(n), bad)
  }
  
  message("Rendering SROC plots for ", K, " diagnostic test(s)...")
  
  # --- graphics: save & protect par; ensure usable device/margins ---
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(cex.axis = cex.axis, cex.lab = cex.lab)
  if (K <= 2) {
    par(mar = c(4, 4, 2, 1) + 0.1, mgp = c(2, 0.6, 0), tcl = -0.3)
  } else if (K == 3) {
    par(mar = c(3.6, 3.6, 1.6, 0.8) + 0.1, mgp = c(2, 0.6, 0), tcl = -0.3)
  } else {
    par(mar = c(3.2, 3.2, 1.2, 0.6) + 0.1, mgp = c(2, 0.6, 0), tcl = -0.3)
  }

  
  if (K == 1) {
    par(mfrow = c(1, 2))
  } else if (K == 2) {
    par(mfrow = c(1, 3))
  } else if (K == 3) {
    par(mfrow = c(2, 2))
  } else if (K >= 4) {
    par(mfrow = c(3, 2))
  }
  
  ## ---------------- K == 1 ----------------
  if (K == 1) {
    sp.range <- matrix(range(stud.sp[, 1]), ncol = K, nrow = 2, byrow = FALSE)
    pool.se  <- result["post.Se[1]", 3]
    pool.sp  <- result["post.Sp[1]", 3]
    alpha.post <- rbind(samp[[1]][, "mu[2]"])
    beta.post  <- rbind(samp[[1]][, "mu[3]"])
    
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    for (i in seq_along(x1)) {
      y_1l[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[6]],  cov[[9]],  0.025)
      y_1m[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[6]],  cov[[9]],  0.50)
      y_1u[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[6]],  cov[[9]],  0.975)
    }
    
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", main = testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)),
            density = c(50, 50), angle = c(-45, 45), col = "lightgrey")
    
    legend1 <- c(paste("SROC for ", testname[1]),
                 paste("Estimated points of ", testname[1]),
                 paste("Pooled TP and FP for ", testname[1]))
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", ...)
    keep1 <- safe_keep(l[[1]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch = 1, col = "black")
    points(1 - pool.sp, pool.se, pch = 3, col = "black")
    legend("bottomright", legend = legend1,
           lty = c(1, NA, NA), pch = c(NA, 1, 3),
           col = c("black", "black", "black"), cex = 0.5)
  }
  
  ## ---------------- K == 2 ----------------
  if (K == 2) {
    sp.range <- matrix(c(range(stud.sp[, 1]), range(stud.sp[, 2])),
                       ncol = K, nrow = 2, byrow = FALSE)
    pool.se <- result[c("post.Se[1]", "post.Se[2]"), 3]
    pool.sp <- result[c("post.Sp[1]", "post.Sp[2]"), 3]
    alpha.post <- rbind(samp[[1]][, "mu[2]"], samp[[1]][, "mu[4]"])
    beta.post  <- rbind(samp[[1]][, "mu[3]"], samp[[1]][, "mu[5]"])
    
    # test 1
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    for (i in seq_along(x1)) {
      y_1l[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[8]],  cov[[13]], 0.025)
      y_1m[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[8]],  cov[[13]], 0.50)
      y_1u[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[8]],  cov[[13]], 0.975)
    }
    # test 2
    x2 <- seq(1 - sp.range[2, 2], 1 - sp.range[1, 2], length.out = 2000)
    y_2l <- y_2m <- y_2u <- numeric(length(x2))
    for (i in seq_along(x2)) {
      y_2l[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[20]], cov[[25]], 0.025)
      y_2m[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[20]], cov[[25]], 0.50)
      y_2u[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[20]], cov[[25]], 0.975)
    }
    
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", main = testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)), 
            density = c(50, 50), angle = c(-45, 45),
            col = "lightgrey", border = NA)
    
    plot(x2, y_2m, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "blue",
         xlab = "False positive rate", ylab = "True positive rate", main = testname[2], ...)
    polygon(c(x2, rev(x2)), c(y_2l, rev(y_2u)), 
            density = c(50, 50), angle = c(-45, 45),
            col = "lightgrey", border = NA)
    
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", ...)
    lines(x2, y_2m, lty = 2, col = "blue")
    
    legend2 <- c(paste("SROC for ", testname[1]),
                 paste("SROC for ", testname[2]),
                 paste("Estimated points of ", testname[1]),
                 paste("Estimated points of ", testname[2]),
                 paste("Pooled TP and FP for ", testname[1]),
                 paste("Pooled TP and FP for ", testname[2]))
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch = 1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch = 2, col = "blue")
    points(1 - pool.sp, pool.se, pch = c(3, 4), col = c("black", "blue"))
    legend("bottomright", legend = legend2,
           lty = c(1, 2, NA, NA, NA, NA), pch = c(NA, NA, 1, 2, 3, 4),
           col = c("black", "blue", "black", "blue", "black", "blue"),
           cex = 0.55)
  }
  
  ## ---------------- K == 3 ----------------
  if (K == 3) {
    sp.range <- matrix(c(range(stud.sp[, 1]), range(stud.sp[, 2]), range(stud.sp[, 3])),
                       ncol = K, nrow = 2, byrow = FALSE)
    pool.se <- result[c("post.Se[1]","post.Se[2]","post.Se[3]"), 3]
    pool.sp <- result[c("post.Sp[1]","post.Sp[2]","post.Sp[3]"), 3]
    alpha.post <- rbind(samp[[1]][, "mu[2]"], samp[[1]][, "mu[4]"], samp[[1]][, "mu[6]"])
    beta.post  <- rbind(samp[[1]][, "mu[3]"], samp[[1]][, "mu[5]"], samp[[1]][, "mu[7]"])
    
    # test 1
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    for (i in seq_along(x1)) {
      y_1l[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[10]], cov[[17]], 0.025)
      y_1m[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[10]], cov[[17]], 0.50)
      y_1u[i] <- cb_h(x1[i], alpha.post[1, ], beta.post[1, ], cov[[10]], cov[[17]], 0.975)
    }
    # test 2
    x2 <- seq(1 - sp.range[2, 2], 1 - sp.range[1, 2], length.out = 2000)
    y_2l <- y_2m <- y_2u <- numeric(length(x2))
    for (i in seq_along(x2)) {
      y_2l[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[26]], cov[[33]], 0.025)
      y_2m[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[26]], cov[[33]], 0.50)
      y_2u[i] <- cb_h(x2[i], alpha.post[2, ], beta.post[2, ], cov[[26]], cov[[33]], 0.975)
    }
    # test 3
    x3 <- seq(1 - sp.range[2, 3], 1 - sp.range[1, 3], length.out = 2000)
    y_3l <- y_3m <- y_3u <- numeric(length(x3))
    for (i in seq_along(x3)) {
      y_3l[i] <- cb_h(x3[i], alpha.post[3, ], beta.post[3, ], cov[[42]], cov[[49]], 0.025)
      y_3m[i] <- cb_h(x3[i], alpha.post[3, ], beta.post[3, ], cov[[42]], cov[[49]], 0.50)
      y_3u[i] <- cb_h(x3[i], alpha.post[3, ], beta.post[3, ], cov[[42]], cov[[49]], 0.975)
    }
    
    par(mfrow = c(2, 2))
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", main = testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)),
            density = c(50, 50), angle = c(-45, 45), col = "lightgrey")
    
    plot(x2, y_2m, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "blue",
         xlab = "False positive rate", ylab = "True positive rate", main = testname[2], ...)
    polygon(c(x2, rev(x2)), c(y_2l, rev(y_2u)),
            density = c(50, 50), angle = c(-45, 45), col = "lightgrey")
    
    plot(x3, y_3m, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "blue",
         xlab = "False positive rate", ylab = "True positive rate", main = testname[3], ...)
    polygon(c(x3, rev(x3)), c(y_3l, rev(y_3u)),
            density = c(50, 50), angle = c(-45, 45), col = "lightgrey")
    
    plot(x1, y_1m, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False positive rate", ylab = "True positive rate", ...)
    lines(x2, y_2m, lty = 2, col = "blue")
    lines(x3, y_3m, lty = 3, col = "green")
    
    legend3 <- c(paste("SROC for ", testname[1]),
                 paste("SROC for ", testname[2]),
                 paste("SROC for ", testname[3]),
                 paste("Estimated points of ", testname[1]),
                 paste("Estimated points of ", testname[2]),
                 paste("Estimated points of ", testname[3]),
                 paste("Pooled TP and FP for ", testname[1]),
                 paste("Pooled TP and FP for ", testname[2]),
                 paste("Pooled TP and FP for ", testname[3]))
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    keep3 <- safe_keep(l[[3]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch = 1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch = 2, col = "blue")
    points(1 - stud.sp[keep3, 3], stud.se[keep3, 3], pch = 3, col = "green")
    points(1 - pool.sp, pool.se, pch = c(4, 5, 6), col = c("black", "blue", "green"))
    legend("bottomright", legend = legend3,
           lty = c(1, 2, 3, NA, NA, NA, NA, NA, NA),
           pch = c(NA, NA, NA, 1, 2, 3, 4, 5, 6),
           col = c("black", "blue", "green", "black", "blue", "green", "black", "blue", "green"),
           cex = 0.35)
  }
  
  ## ---------------- K == 4 ----------------
  if (K == 4) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2]),
                         range(stud.sp[,3]), range(stud.sp[,4])),
                       ncol = K, nrow = 2, byrow = FALSE)
    alpha.post <- rbind(samp[[1]][, "mu[2]"], samp[[1]][, "mu[4]"],
                        samp[[1]][, "mu[6]"], samp[[1]][, "mu[8]"])
    beta.post  <- rbind(samp[[1]][, "mu[3]"], samp[[1]][, "mu[5]"],
                        samp[[1]][, "mu[7]"], samp[[1]][, "mu[9]"])
    pool.se <- result[c("post.Se[1]","post.Se[2]","post.Se[3]","post.Se[4]"), 3]
    pool.sp <- result[c("post.Sp[1]","post.Sp[2]","post.Sp[3]","post.Sp[4]"), 3]
    
    # x/y for 1..4
    build_xy <- function(j, covL, covU) {
      x  <- seq(1 - sp.range[2, j], 1 - sp.range[1, j], length.out = 2000)
      yl <- ym <- yu <- numeric(length(x))
      for (i in seq_along(x)) {
        yl[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.025)
        ym[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.50)
        yu[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.975)
      }
      list(x=x, yl=yl, ym=ym, yu=yu)
    }
    xy1 <- build_xy(1, 12, 21)
    xy2 <- build_xy(2, 32, 41)
    xy3 <- build_xy(3, 52, 61)
    xy4 <- build_xy(4, 72, 81)
    
    par(mfrow = c(3, 2))
    plot(xy1$x, xy1$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", main=testname[1], ...)
    polygon(c(xy1$x, rev(xy1$x)), c(xy1$yl, rev(xy1$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy2$x, xy2$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[2], ...)
    polygon(c(xy2$x, rev(xy2$x)), c(xy2$yl, rev(xy2$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy3$x, xy3$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[3], ...)
    polygon(c(xy3$x, rev(xy3$x)), c(xy3$yl, rev(xy3$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy4$x, xy4$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[4], ...)
    polygon(c(xy4$x, rev(xy4$x)), c(xy4$yl, rev(xy4$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy1$x, xy1$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(xy2$x, xy2$ym, lty=2, col="blue")
    lines(xy3$x, xy3$ym, lty=3, col="green")
    lines(xy4$x, xy4$ym, lty=4, col="red")
    
    legend4 <- c(paste("SROC for ", testname[1]),
                 paste("SROC for ", testname[2]),
                 paste("SROC for ", testname[3]),
                 paste("SROC for ", testname[4]),
                 paste("Estimated points of ", testname[1]),
                 paste("Estimated points of ", testname[2]),
                 paste("Estimated points of ", testname[3]),
                 paste("Estimated points of ", testname[4]),
                 paste("Pooled TP and FP for ", testname[1]),
                 paste("Pooled TP and FP for ", testname[2]),
                 paste("Pooled TP and FP for ", testname[3]),
                 paste("Pooled TP and FP for ", testname[4]))
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    keep3 <- safe_keep(l[[3]], nstu)
    keep4 <- safe_keep(l[[4]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch=1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch=2, col="blue")
    points(1 - stud.sp[keep3, 3], stud.se[keep3, 3], pch=3, col="green")
    points(1 - stud.sp[keep4, 4], stud.se[keep4, 4], pch=4, col="red")
    points(1 - pool.sp, pool.se, pch=c(5,6,7,8), col=c("black","blue","green","red"))
    legend("bottomright", legend=legend4,
           lty=c(1,2,3,4,NA,NA,NA,NA,NA,NA,NA,NA),
           pch=c(NA,NA,NA,NA,1,2,3,4,5,6,7,8),
           col=c("black","blue","green","red","black","blue","green","red","black","blue","green","red"),
           cex=0.35)
  }
  
  ## ---------------- K == 5 ----------------
  if (K == 5) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2]),
                         range(stud.sp[,3]), range(stud.sp[,4]), range(stud.sp[,5])),
                       ncol = K, nrow = 2, byrow = FALSE)
    alpha.post <- rbind(samp[[1]][, "mu[2]"],  samp[[1]][, "mu[4]"],
                        samp[[1]][, "mu[6]"],  samp[[1]][, "mu[8]"],
                        samp[[1]][, "mu[10]"])
    beta.post  <- rbind(samp[[1]][, "mu[3]"],  samp[[1]][, "mu[5]"],
                        samp[[1]][, "mu[7]"],  samp[[1]][, "mu[9]"],
                        samp[[1]][, "mu[11]"])
    pool.se <- result[paste0("post.Se[", 1:5, "]"), 3]
    pool.sp <- result[paste0("post.Sp[", 1:5, "]"), 3]
    
    build_xy <- function(j, covL, covU) {
      x  <- seq(1 - sp.range[2, j], 1 - sp.range[1, j], length.out = 2000)
      yl <- ym <- yu <- numeric(length(x))
      for (i in seq_along(x)) {
        yl[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.025)
        ym[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.50)
        yu[i] <- cb_h(x[i], alpha.post[j, ], beta.post[j, ], cov[[covL]], cov[[covU]], 0.975)
      }
      list(x=x, yl=yl, ym=ym, yu=yu)
    }
    xy1 <- build_xy(1, 14, 25)
    xy2 <- build_xy(2, 38, 49)
    xy3 <- build_xy(3, 62, 73)
    xy4 <- build_xy(4, 86, 97)
    xy5 <- build_xy(5, 110, 121)
    
    par(mfrow = c(3, 2))
    plot(xy1$x, xy1$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", main=testname[1], ...)
    polygon(c(xy1$x, rev(xy1$x)), c(xy1$yl, rev(xy1$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy2$x, xy2$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[2], ...)
    polygon(c(xy2$x, rev(xy2$x)), c(xy2$yl, rev(xy2$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy3$x, xy3$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[3], ...)
    polygon(c(xy3$x, rev(xy3$x)), c(xy3$yl, rev(xy3$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy4$x, xy4$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[4], ...)
    polygon(c(xy4$x, rev(xy4$x)), c(xy4$yl, rev(xy4$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy5$x, xy5$ym, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[5], ...)
    polygon(c(xy5$x, rev(xy5$x)), c(xy5$yl, rev(xy5$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(xy1$x, xy1$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(xy2$x, xy2$ym, lty=2, col="blue")
    lines(xy3$x, xy3$ym, lty=3, col="green")
    lines(xy4$x, xy4$ym, lty=4, col="red")
    lines(xy5$x, xy5$ym, lty=5, col="yellow")
    
    legend5 <- c(paste("SROC for ", testname[1]),
                 paste("SROC for ", testname[2]),
                 paste("SROC for ", testname[3]),
                 paste("SROC for ", testname[4]),
                 paste("SROC for ", testname[5]),
                 paste("Estimated points of ", testname[1]),
                 paste("Estimated points of ", testname[2]),
                 paste("Estimated points of ", testname[3]),
                 paste("Estimated points of ", testname[4]),
                 paste("Estimated points of ", testname[5]),
                 paste("Pooled TP and FP for ", testname[1]),
                 paste("Pooled TP and FP for ", testname[2]),
                 paste("Pooled TP and FP for ", testname[3]),
                 paste("Pooled TP and FP for ", testname[4]),
                 paste("Pooled TP and FP for ", testname[5]))
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    keep3 <- safe_keep(l[[3]], nstu)
    keep4 <- safe_keep(l[[4]], nstu)
    keep5 <- safe_keep(l[[5]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch=1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch=2, col="blue")
    points(1 - stud.sp[keep3, 3], stud.se[keep3, 3], pch=3, col="green")
    points(1 - stud.sp[keep4, 4], stud.se[keep4, 4], pch=4, col="red")
    points(1 - stud.sp[keep5, 5], stud.se[keep5, 5], pch=5, col="yellow")
    points(1 - pool.sp, pool.se, pch=c(6,7,8,9,10), col=c("black","blue","green","red","yellow"))
    legend("bottomright", legend=legend5,
           lty=c(1,2,3,4,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
           pch=c(NA,NA,NA,NA,NA,1,2,3,4,5,6,7,8,9,10),
           col=c("black","blue","green","red","yellow",
                 "black","blue","green","red","yellow",
                 "black","blue","green","red","yellow"),
           cex=0.35)
  }
  
  invisible(NULL)
}
