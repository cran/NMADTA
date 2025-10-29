
sroc.hsroc <- function(x,
                       cex.axis = 1.2,
                       cex.lab  = 1.2,
                       ...) {
  # ---- Checks ----
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hsroc")
    stop("sroc.hsroc() is for HSROC models only.")
  samp <- x$rawOutput
  dat  <- x$dat
  nstu <- x$nstu
  K    <- x$K
  if (!is.list(samp))
    stop("x$rawOutput must contain a list of MCMC chains (list).")
  
  # ---- test names ----
  if (!is.null(x$testname)) {
    testname <- x$testname
  } else if (!is.null(x$testnames)) {
    testname <- x$testnames
  } else {
    testname <- paste0("Test ", seq_len(K))
  }
  if (length(testname) != K) testname <- paste0("Test ", seq_len(K))
  
  # ---- Core prep ----
  index <- indicator(K, nstu, dat)
  l      <- except(index, K, nstu)
  summ   <- summary(samp)
  result <- summ[[2]]
  
  # HSROC
  para_study_se <- paralist("stud.se", summ)
  para_study_sp <- paralist("stud.sp", summ)
  stud.se <- matrix(result[para_study_se, 3], ncol = K, nrow = nstu, byrow = FALSE)
  stud.sp <- matrix(result[para_study_sp, 3], ncol = K, nrow = nstu, byrow = FALSE)
  
  # pooled
  pool.se <- result[paste0("post.se[", seq_len(K), "]"), 3]
  pool.sp <- result[paste0("post.sp[", seq_len(K), "]"), 3]
  
  ###
  pull_chain <- function(par_name) {
    unlist(lapply(samp, function(ch) ch[, par_name]), use.names = FALSE)
  }
  beta.post <- do.call(rbind, lapply(seq_len(K), function(k) pull_chain(sprintf("beta[%d]", k))))
  rownames(beta.post) <- paste0("beta", seq_len(K))
  mu2.post  <- do.call(rbind, lapply(seq_len(K), function(k) pull_chain(sprintf("mu2[%d]", k))))
  rownames(mu2.post)  <- paste0("mu2",  seq_len(K))
  
  # helpers
  safe_keep <- function(lj, n) {
    bad <- unlist(lj, use.names = FALSE)
    if (length(bad) == 0) return(seq_len(n))
    setdiff(seq_len(n), bad)
  }
  
  message("Rendering HSROC SROC plots for ", K, " diagnostic test(s)...")
  
  # ---- graphics: device & par ----
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
  
  ## ------------- K == 1 -------------
  if (K == 1) {
    sp.range <- matrix(range(stud.sp[, 1]), ncol = K, nrow = 2, byrow = FALSE)
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    for (i in seq_along(x1)) {
      y_1l[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.025)
      y_1m[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.50)
      y_1u[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.975)
    }
    
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate",
         main=testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    keep1 <- safe_keep(l[[1]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch=1, col="black")
    points(1 - pool.sp[1], pool.se[1], pch=3, col="black")
    legend("bottomright",
           legend=c(paste("SROC for ", testname[1]),
                    paste("Estimated points of ", testname[1]),
                    paste("Pooled TP and FP for ", testname[1])),
           lty=c(1, NA, NA), pch=c(NA, 1, 3), col="black", cex=0.8)
  }
  
  ## ------------- K == 2 -------------
  if (K == 2) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2])),
                       ncol = K, nrow = 2, byrow = FALSE)
    
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    x2 <- seq(1 - sp.range[2, 2], 1 - sp.range[1, 2], length.out = 2000)
    
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    y_2l <- y_2m <- y_2u <- numeric(length(x2))
    
    for (i in seq_along(x1)) {
      y_1l[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.025)
      y_1m[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.50)
      y_1u[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.975)
    }
    for (i in seq_along(x2)) {        
      y_2l[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.025)
      y_2m[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.50)
      y_2u[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.975)
    }
    
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate",
         main=testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)),
            density=c(50,50), angle=c(-45,45),
            col="lightgrey", border=NA)
    
    plot(x2, y_2m, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate",
         main=testname[2], ...)
    polygon(c(x2, rev(x2)), c(y_2l, rev(y_2u)),
            density=c(50,50), angle=c(-45,45),
            col="lightgrey", border=NA)
    
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(x2, y_2m, lty=2, col="blue")
    
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch=1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch=2, col="blue")
    points(1 - pool.sp, pool.se, pch=c(3,4), col=c("black","blue"))
    legend("bottomright",
           legend=c(paste("SROC for ", testname[1]),
                    paste("SROC for ", testname[2]),
                    paste("Estimated points of ", testname[1]),
                    paste("Estimated points of ", testname[2]),
                    paste("Pooled TP and FP for ", testname[1]),
                    paste("Pooled TP and FP for ", testname[2])),
           lty=c(1,2,NA,NA,NA,NA),
           pch=c(NA,NA,1,2,3,4),
           col=c("black","blue","black","blue","black","blue"),
           cex=0.6)
  }
  
  ## ------------- K == 3 -------------
  if (K == 3) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2]), range(stud.sp[,3])),
                       ncol = K, nrow = 2, byrow = FALSE)
    x1 <- seq(1 - sp.range[2, 1], 1 - sp.range[1, 1], length.out = 2000)
    x2 <- seq(1 - sp.range[2, 2], 1 - sp.range[1, 2], length.out = 2000)
    x3 <- seq(1 - sp.range[2, 3], 1 - sp.range[1, 3], length.out = 2000)
    
    y_1l <- y_1m <- y_1u <- numeric(length(x1))
    y_2l <- y_2m <- y_2u <- numeric(length(x2))
    y_3l <- y_3m <- y_3u <- numeric(length(x3))
    
    for (i in seq_along(x1)) {
      y_1l[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.025)
      y_1m[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.50)
      y_1u[i] <- cb(x1[i], beta.post[1, ], mu2.post[1, ], 0.975)
    }
    for (i in seq_along(x2)) {
      y_2l[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.025)
      y_2m[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.50)
      y_2u[i] <- cb(x2[i], beta.post[2, ], mu2.post[2, ], 0.975)
    }
    for (i in seq_along(x3)) {
      y_3l[i] <- cb(x3[i], beta.post[3, ], mu2.post[3, ], 0.025)
      y_3m[i] <- cb(x3[i], beta.post[3, ], mu2.post[3, ], 0.50)
      y_3u[i] <- cb(x3[i], beta.post[3, ], mu2.post[3, ], 0.975)
    }
    
    par(mfrow = c(2, 2))
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", main=testname[1], ...)
    polygon(c(x1, rev(x1)), c(y_1l, rev(y_1u)), col="lightgrey", border=NA)
    
    plot(x2, y_2m, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[2], ...)
    polygon(c(x2, rev(x2)), c(y_2l, rev(y_2u)), col="lightgrey", border=NA)
    
    plot(x3, y_3m, type="l", xlim=c(0,1), ylim=c(0,1), col="blue",
         xlab="False positive rate", ylab="True positive rate", main=testname[3], ...)
    polygon(c(x3, rev(x3)), c(y_3l, rev(y_3u)), col="lightgrey", border=NA)
    
    plot(x1, y_1m, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(x2, y_2m, lty=2, col="blue")
    lines(x3, y_3m, lty=3, col="green")
    
    keep1 <- safe_keep(l[[1]], nstu)
    keep2 <- safe_keep(l[[2]], nstu)
    keep3 <- safe_keep(l[[3]], nstu)
    points(1 - stud.sp[keep1, 1], stud.se[keep1, 1], pch=1)
    points(1 - stud.sp[keep2, 2], stud.se[keep2, 2], pch=2, col="blue")
    points(1 - stud.sp[keep3, 3], stud.se[keep3, 3], pch=3, col="green")
    points(1 - pool.sp, pool.se, pch=c(4,5,6), col=c("black","blue","green"))
    legend("bottomright",
           legend=c(paste("SROC for ", testname[1]),
                    paste("SROC for ", testname[2]),
                    paste("SROC for ", testname[3]),
                    paste("Estimated points of ", testname[1]),
                    paste("Estimated points of ", testname[2]),
                    paste("Estimated points of ", testname[3]),
                    paste("Pooled TP and FP for ", testname[1]),
                    paste("Pooled TP and FP for ", testname[2]),
                    paste("Pooled TP and FP for ", testname[3])),
           lty=c(1,2,3,NA,NA,NA,NA,NA,NA),
           pch=c(NA,NA,NA,1,2,3,4,5,6),
           col=c("black","blue","green","black","blue","green","black","blue","green"),
           cex=0.45)
  }
  
  ## ------------- K == 4 -------------
  if (K == 4) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2]),
                         range(stud.sp[,3]), range(stud.sp[,4])),
                       ncol = K, nrow = 2, byrow = FALSE)
    
    build_xy <- function(j) {
      x  <- seq(1 - sp.range[2, j], 1 - sp.range[1, j], length.out = 2000)
      yl <- ym <- yu <- numeric(length(x))
      for (i in seq_along(x)) {
        yl[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.025)
        ym[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.50)
        yu[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.975)
      }
      list(x=x, yl=yl, ym=ym, yu=yu)
    }
    xy <- lapply(1:4, build_xy)
    
    par(mfrow = c(3, 2))
    cols <- c("black", "blue", "blue", "blue")
    for (j in 1:4) {
      plot(xy[[j]]$x, xy[[j]]$ym, type="l", xlim=c(0,1), ylim=c(0,1),
           xlab="False positive rate", ylab="True positive rate",
           main=testname[j], col=cols[j], ...)
      polygon(c(xy[[j]]$x, rev(xy[[j]]$x)), c(xy[[j]]$yl, rev(xy[[j]]$yu)),
              density=c(50,50), angle=c(-45,45), col="lightgrey")
      if (j == 2 || j == 3 || j == 4) {} # 仅用于对齐结构
    }
    
    # 综合图
    plot(xy[[1]]$x, xy[[1]]$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(xy[[2]]$x, xy[[2]]$ym, lty=2, col="blue")
    lines(xy[[3]]$x, xy[[3]]$ym, lty=3, col="green")
    lines(xy[[4]]$x, xy[[4]]$ym, lty=4, col="red")
    
    keep <- lapply(1:4, function(j) safe_keep(l[[j]], nstu))
    points(1 - stud.sp[keep[[1]], 1], stud.se[keep[[1]], 1], pch=1)
    points(1 - stud.sp[keep[[2]], 2], stud.se[keep[[2]], 2], pch=2, col="blue")
    points(1 - stud.sp[keep[[3]], 3], stud.se[keep[[3]], 3], pch=3, col="green")
    points(1 - stud.sp[keep[[4]], 4], stud.se[keep[[4]], 4], pch=4, col="red")
    points(1 - pool.sp, pool.se, pch=c(5,6,7,8), col=c("black","blue","green","red"))
    legend("bottomright",
           legend=c(paste("SROC for ", testname[1]),
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
                    paste("Pooled TP and FP for ", testname[4])),
           lty=c(1,2,3,4,NA,NA,NA,NA,NA,NA,NA,NA),
           pch=c(NA,NA,NA,NA,1,2,3,4,5,6,7,8),
           col=c("black","blue","green","red","black","blue","green","red","black","blue","green","red"),
           cex=0.35)
  }
  
  ## ------------- K == 5 -------------
  if (K == 5) {
    sp.range <- matrix(c(range(stud.sp[,1]), range(stud.sp[,2]),
                         range(stud.sp[,3]), range(stud.sp[,4]),
                         range(stud.sp[,5])),
                       ncol = K, nrow = 2, byrow = FALSE)
    
    build_xy <- function(j) {
      x  <- seq(1 - sp.range[2, j], 1 - sp.range[1, j], length.out = 2000)
      yl <- ym <- yu <- numeric(length(x))
      for (i in seq_along(x)) {
        yl[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.025)
        ym[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.50)
        yu[i] <- cb(x[i], beta.post[j, ], mu2.post[j, ], 0.975)
      }
      list(x=x, yl=yl, ym=ym, yu=yu)
    }
    xy <- lapply(1:5, build_xy)
    
    par(mfrow = c(3, 2))
    cols <- c("black", "blue", "blue", "blue", "blue")
    for (j in 1:5) {
      if (j <= 5 && j <= 4) {
        plot(xy[[j]]$x, xy[[j]]$ym, type="l", xlim=c(0,1), ylim=c(0,1),
             xlab="False positive rate", ylab="True positive rate",
             main=testname[j], col=cols[j], ...)
        polygon(c(xy[[j]]$x, rev(xy[[j]]$x)), c(xy[[j]]$yl, rev(xy[[j]]$yu)),
                density=c(50,50), angle=c(-45,45), col="lightgrey")
      }
    }
    # 第5张单图
    plot(xy[[5]]$x, xy[[5]]$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate",
         main=testname[5], col="blue", ...)
    polygon(c(xy[[5]]$x, rev(xy[[5]]$x)), c(xy[[5]]$yl, rev(xy[[5]]$yu)),
            density=c(50,50), angle=c(-45,45), col="lightgrey")
    
    # 综合图
    plot(xy[[1]]$x, xy[[1]]$ym, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate", ylab="True positive rate", ...)
    lines(xy[[2]]$x, xy[[2]]$ym, lty=2, col="blue")
    lines(xy[[3]]$x, xy[[3]]$ym, lty=3, col="green")
    lines(xy[[4]]$x, xy[[4]]$ym, lty=4, col="red")
    lines(xy[[5]]$x, xy[[5]]$ym, lty=5, col="yellow")
    
    keep <- lapply(1:5, function(j) safe_keep(l[[j]], nstu))
    points(1 - stud.sp[keep[[1]], 1], stud.se[keep[[1]], 1], pch=1)
    points(1 - stud.sp[keep[[2]], 2], stud.se[keep[[2]], 2], pch=2, col="blue")
    points(1 - stud.sp[keep[[3]], 3], stud.se[keep[[3]], 3], pch=3, col="green")
    points(1 - stud.sp[keep[[4]], 4], stud.se[keep[[4]], 4], pch=4, col="red")
    points(1 - stud.sp[keep[[5]], 5], stud.se[keep[[5]], 5], pch=5, col="yellow")
    points(1 - pool.sp, pool.se, pch=c(6,7,8,9,10), col=c("black","blue","green","red","yellow"))
    legend("bottomright",
           legend=c(paste("SROC for ", testname[1]),
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
                    paste("Pooled TP and FP for ", testname[5])),
           lty=c(1,2,3,4,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
           pch=c(NA,NA,NA,NA,NA,1,2,3,4,5,6,7,8,9,10),
           col=c("black","blue","green","red","yellow",
                 "black","blue","green","red","yellow",
                 "black","blue","green","red","yellow"),
           cex=0.35)
  }
  
  invisible(NULL)
}
