#' Summary method for `nmadt` objects
#'
#' Provides a concise summary of posterior results from network meta-analysis
#' of diagnostic tests fitted using the \code{NMADTA} framework.  
#' Displays median estimates and 95% credible intervals for
#' sensitivity, specificity, predictive values, likelihood ratios, and prevalence.
#' 
#' @param object An object of class \code{nmadt}, 
#' typically created by \code{nmadt.hierarchical()} or \code{nmadt.hsroc()}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The function extracts and prints key posterior summaries from the fitted model, 
#' including medians and 95% equal-tail credible intervals (CrI) for diagnostic accuracy parameters.
#' 
#' The output is formatted for human readability, with each section clearly labeled.
#' 
#' @return
#' The function returns the input \code{object} (invisibly) after printing the summary.
#'
#' @examples
#' \donttest{
#' data(dat.kang)
#' set.seed(9)
#' kang.out <- nmadt.hierarchical(nstu=12, K=2, data=dat.kang, 
#'             testname=c("D-dimer","Ultrasonography"))
#' summary(kang.out) 
#' }
#' 
#' @exportS3Method summary nmadt
summary.nmadt <- function(object, ...) {
  cat("Model:", object$model, "\n\n")
  
  cat("Sensitivity (Median and 95% CrI):\n")
  print(object$Se$Median_CI)
  
  cat("\nSpecificity (Median and 95% CrI):\n")
  print(object$Sp$Median_CI)
  
  cat("\nPrevalence:\n")
  print(object$prevalence$Median_CI)
  
  cat("\nPPV:\n")
  print(object$ppv$Median_CI)
  
  cat("\nNPV:\n")
  print(object$npv$Median_CI)
  
  cat("\nLR+:\n")
  print(object$LRpos$Median_CI)
  
  cat("\nLR-:\n")
  print(object$LRneg$Median_CI)
  
  invisible(object)
}