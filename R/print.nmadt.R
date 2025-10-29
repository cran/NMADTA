#' Print method for `nmadt` objects
#'
#' @param x An object of class `nmadt`
#' @param ... Not used. Included for S3 method compatibility.
#'
#' @examples
#' \donttest{
#' data(dat.kang)
#' set.seed(9)
#' kang.out <- nmadt.hierarchical(nstu=12, K=2, data=dat.kang, 
#'             testname=c("D-dimer","Ultrasonography"))
#' print(kang.out) 
#' }
#' 
#' @exportS3Method print nmadt
print.nmadt <- function(x, ...) {
  if (!inherits(x, "nmadt")) {
    stop("Expected an object of class 'nmadt'.")
  }
  
  plusminus <- "\u00B1" 
  
  cat("Model type:", if (!is.null(x$model)) x$model else "NA", "\n\n")
  
  show_blocks <- function(title, comp) {
    cat(sprintf("%s (Mean %s SD):\n", title, plusminus))
    if (!is.null(comp) && !is.null(comp$Mean_SD)) {
      print(comp$Mean_SD)
    } else {
      cat("[not available]\n")
    }
    
    cat(sprintf("\n%s (Median and 95%% CrI):\n", title))
    if (!is.null(comp) && !is.null(comp$Median_CI)) {
      mc <- comp$Median_CI
      cn <- colnames(mc)
      if (!is.null(cn)) {
        colnames(mc) <- gsub("95% CI", "95% CrI", cn, fixed = TRUE)
      }
      print(mc)
    } else {
      cat("[not available]\n")
    }
    cat("\n")
  }
  
  show_blocks("Sensitivity", x$Se)
  show_blocks("Specificity", x$Sp)
  
  invisible(x)
}