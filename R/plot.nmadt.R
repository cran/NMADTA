#' Plot method for `nmadt` objects
#'
#' This method automatically generates diagnostic meta-analysis plots 
#' based on the fitted `nmadt` object and the specified plot `type`.
#'
#' The available plot types are:
#' 
#' \describe{
#'   \item{\strong{"sroc"} (default)}{
#'     Summary Receiver Operating Characteristic (SROC) curve.
#'     Visualizes the trade-off between sensitivity and specificity 
#'     across studies, along with the hierarchical model fit.
#'   }
#'   \item{\strong{"density"}}{
#'     Posterior density plots for study- and test-level sensitivity 
#'     and specificity parameters.
#'     Useful for checking convergence and posterior uncertainty.
#'   }
#'   \item{\strong{"forest"}}{
#'     Forest plot summarizing point estimates and uncertainty intervals 
#'     for sensitivity and specificity of each test across studies.
#'     Helpful for visualizing study heterogeneity.
#'   }
#'   \item{\strong{"contour"}}{
#'     Contour-enhanced plots showing joint posterior density of 
#'     sensitivity and specificity for each test or study.
#'     Useful for visual comparison of test performance.
#'   }
#' }
#'
#' @param x An object of class `nmadt`, typically produced by one of the 
#'   model-fitting functions:
#'   \itemize{
#'     \item \code{nmadt.hierarchical()} — hierarchical model under MAR assumption;
#'     \item \code{nmadt.hsroc()} — HSROC model under MAR assumption;
#'     \item \code{nmadt.hierarchical.MNAR()} — hierarchical model allowing 
#'     for MNAR (missing not at random) mechanism;
#'     \item \code{nmadt.hsroc.MNAR()} — HSROC model allowing for MNAR mechanism.
#'   }
#'   These functions all return an object of class `nmadt` suitable for plotting.
#' @param type Character string specifying the type of plot to generate.
#'   One of \code{"sroc"}, \code{"density"}, \code{"forest"}, or \code{"contour"}.
#'   Defaults to \code{"sroc"} if not specified.
#' @param ... Additional arguments passed to the underlying plotting functions 
#'   (e.g., graphical parameters such as \code{cex.axis}, \code{cex.lab}, etc.).
#'
#' @details
#' If \code{type} is not specified, the function defaults to \code{"sroc"}.
#' For example, both \code{plot(x)} and \code{plot(x, type = "sroc")} 
#' will produce the SROC plot.
#'
#' @method plot nmadt
#' 
#' @examples
#' \donttest{
#' data(dat.kang)
#' set.seed(9)
#' kang.out <- nmadt.hierarchical(nstu=12, K=2, data=dat.kang, 
#'             testname=c("D-dimer","Ultrasonography"))
#' plot(kang.out, type = "sroc")
#' plot(kang.out, type = "forest")
#' plot(kang.out, type = "contour")
#' plot(kang.out, type = "density")
#' }
#' 
#' @return 
#' Invisibly returns the input `nmadt` object \code{x}.  
#' The function is primarily called for its side effect of generating plots 
#' rather than returning a value.
#' @export
#' @importFrom grDevices dev.size dev.new
plot.nmadt <- function(x, type = c("sroc", "density", "forest", "contour"), ...) {
  type <- match.arg(type)
  model_type <- x$model  # "hierarchical" or "HSROC"
  
  if (type == "sroc") {
    if (model_type == "hierarchical") {
      sroc.hierarchical(x, ...)
    } else if (model_type == "HSROC") {
      sroc.hsroc(x, ...)
    }
  } else if (type == "density") {
    if (model_type == "hierarchical") {
      density_hierarchical(x, ...)
    } else if (model_type == "HSROC") {
      density_hsroc(x, ...)
    }
  } else if (type == "forest") {
    if (model_type == "hierarchical") {
      forestplot.hierarchical(x, ...)
    } else {
      forestplot.hsroc(x, ...)
    }
  } else if (type == "contour") {
    if (model_type == "hierarchical") {
      contour_hierarchical(x, ...)
    } else if (model_type == "HSROC") {
      contour_hsroc(x, ...)
    }
  } else {
    stop("Unknown plot type: choose from 'density', 'forest', 'sroc', or 'contour'")
  }
}