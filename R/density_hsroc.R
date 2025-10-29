
density_hsroc <- function(x,
                          cex.axis = 1.2, cex.lab = 1.2,
                          max_per_page = 3,
                          ...) {
  if (!inherits(x, "nmadt"))
    stop("Input must be an 'nmadt' object.")
  if (tolower(x$model) != "hsroc")
    stop("density_hsroc() is for HSROC models only.")
  
  samp.gen  <- x$rawOutput
  testname  <- x$testname
  K         <- x$K
  n.chains  <- length(samp.gen)
  
  message("Start saving posterior density plots...\n")
  
  # 初始化
  mcmcSe <- vector("list", K)
  mcmcSp <- vector("list", K)
  dens.Se <- matrix(0, K, 3, dimnames = list(NULL, c("ymax","xmin","xmax")))
  dens.Sp <- matrix(0, K, 3, dimnames = list(NULL, c("ymax","xmin","xmax")))
  
  # 自动检测变量命名（兼容 post.se / post.Se）
  var_se_prefix <- if (any(grepl("^post.Se", colnames(samp.gen[[1]])))) "post.Se" else "post.se"
  var_sp_prefix <- if (any(grepl("^post.Sp", colnames(samp.gen[[1]])))) "post.Sp" else "post.sp"
  
  # 提取 MCMC 样本并计算密度范围
  for (i in seq_len(K)) {
    tempSe <- tempSp <- NULL
    for (j in seq_len(n.chains)) {
      tempSe <- c(tempSe, as.vector(samp.gen[[j]][, paste0(var_se_prefix, "[", i, "]")]))
      tempSp <- c(tempSp, as.vector(samp.gen[[j]][, paste0(var_sp_prefix, "[", i, "]")]))
    }
    mcmcSe[[i]] <- tempSe
    mcmcSp[[i]] <- tempSp
    
    tempdens.Se <- density(tempSe)
    dens.Se[i, ] <- c(max(tempdens.Se$y),
                      quantile(tempSe, c(0.001, 0.999)))
    
    tempdens.Sp <- density(tempSp)
    dens.Sp[i, ] <- c(max(tempdens.Sp$y),
                      quantile(tempSp, c(0.001, 0.999)))
  }
  
  ymaxSe <- max(dens.Se[, "ymax"])
  xminSe <- min(dens.Se[, "xmin"]); xmaxSe <- max(dens.Se[, "xmax"])
  ymaxSp <- max(dens.Sp[, "ymax"])
  xminSp <- min(dens.Sp[, "xmin"]); xmaxSp <- max(dens.Sp[, "xmax"])
  
  # 自动分页绘图
  pages <- ceiling(K / max_per_page)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  for (p in seq_len(pages)) {
    idx_start <- (p - 1) * max_per_page + 1
    idx_end   <- min(p * max_per_page, K)
    idx <- idx_start:idx_end
    n_this <- length(idx)
    
    if (n_this > 2) {
      margins <- c(4, 4, 2, 1) + 0.1
      cex.axis <- 0.9; cex.lab <- 0.9
    } else {
      margins <- c(5, 5, 2, 2) + 0.1
    }
    
    par(mfrow = c(n_this, 2), mar = margins)
    
    for (i in idx) {
      plot(density(mcmcSe[[i]]),
           xlim = c(xminSe, xmaxSe), ylim = c(0, ymaxSe),
           xlab = paste(testname[i], " Se(%)"), ylab = "Density",
           main = "", lty = 1, lwd = 2,
           cex.axis = cex.axis, cex.lab = cex.lab)
      
      plot(density(mcmcSp[[i]]),
           xlim = c(xminSp, xmaxSp), ylim = c(0, ymaxSp),
           xlab = paste(testname[i], " Sp(%)"), ylab = "Density",
           main = "", lty = 1, lwd = 2,
           cex.axis = cex.axis, cex.lab = cex.lab)
    }
    
    if (pages > 1 && p < pages) {
      readline(prompt = sprintf("Page %d/%d done. Press [Enter] to continue...", p, pages))
    }
  }

  invisible(recordPlot())  # 返回绘图记录以便保存
}
