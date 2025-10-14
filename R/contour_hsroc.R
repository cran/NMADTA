#' @importFrom ks Hscv kde
#' @import ggplot2
#' @import reshape2
#' @import grDevices
contour_hsroc <- function(samp.gen, K, testname, dirc) {
  plots <- vector("list", K)
  
  grDevices::pdf(file.path(dirc, "contour_hsroc.pdf"))
  
  for (i in seq_len(K)) {
    # Extract posterior sensitivity/specificity
    se <- as.numeric(unlist(lapply(samp.gen, function(x) x[, paste0("post.se[", i, "]")])))
    sp <- as.numeric(unlist(lapply(samp.gen, function(x) x[, paste0("post.sp[", i, "]")])))
    fpr <- 1 - sp
    dat <- data.frame(x = fpr, y = se)
    
    # Kernel density estimation
    hscv <- ks::Hscv(dat)
    fhat <- ks::kde(dat, H = hscv, compute.cont = TRUE)
    dimnames(fhat$estimate) <- list(fhat$eval.points[[1]], fhat$eval.points[[2]])
    
    # Reshape density estimates
    aa <- reshape2::melt(fhat$estimate)
    aa$Var1 <- rep(fhat$eval.points[[1]], times = length(fhat$eval.points[[2]]))
    aa$Var2 <- rep(fhat$eval.points[[2]], each  = length(fhat$eval.points[[1]]))
    
    # Build contour plot
    plots[[i]] <- ggplot(aa, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = value)) +
      geom_contour(aes(z = value), breaks = fhat$cont["75%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["50%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["25%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["10%"]) +
      geom_contour(aes(z = value), breaks = fhat$cont["5%"]) +
      labs(title = testname[i], x = "False positive rate", y = "True positive rate")
    
    print(plots[[i]])
  }
  
  grDevices::dev.off()
  invisible(plots)
}