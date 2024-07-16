#' @title Quasi-likelihood under the Independence Model Criterion
#'
#' @description Calculating the quasi-likelihood under the independence model
#' criterion (QIC) based on Pan (2001a).
#'
#' @param object An object of class \code{"gee"} representing the fit.
#' @param type A character string specifying the type of penalty to calculate.
#' The following are permitted: \code{"pan"}, \code{"hh"}, \code{"u"}.
#' Default is \code{"pan"}.
#' @param tol A numeric value specifying the tolerance when calculating the
#' inverse of the model-based covariance matrix.
#' Default is \code{.Machine$double.eps}.
#' @param envir An environment object to which find the \code{"gee"} object.
#'
#' @details The calculation of QICu proposed by Pan(2001) is calculated when it
#' is assumed that they working correlation structure is properly specified.
#' This value can be obtained by specifying \code{type = "u"}. Specifying
#' \code{type = "hh"} calculates the penalty based on Hardin and Hilbe (2012).
#' More information on this penalty can be found at \code{?gee::cic}.
#'
#' @return A named numeric vector containing the QIC value for the \code{"gee"}
#' object under consideration. The name of the value is based on \code{type},
#' where \code{"QIC"} denotes \code{type = "pan"}, \code{"QIChh"} denotes
#' \code{type = "hh"}, and \code{"QICu"} denotes \code{type = "u"}.
#'
#' @seealso \code{\link{gee::cic}}, \code{\link{gee::quasilik}}
#'
#' @references
#' Pan, W. (2001a). \emph{Akaike's information criterion in generalized
#' estimating equations}. Biometrics, 57, 120-125.\cr
#' Hardin, J.W.  and Hilbe, J.M. (2012). \emph{Generalized
#' Estimating Equations, 2nd Edition}, Chapman and Hall/CRC: New York. \cr
#' Hin, L.-Y. and Wang, Y-G. (2009). \emph{Working-correlation-structure
#' identification in generalized estimating equations},
#' Statistics in Medicine 28: 642-658.
#'
#' @examples
#' # Marginal analysis of random effects model for wool
#' data(warpbreaks)
#' fit <- gee::gee(
#'  breaks ~ tension, id = wool, data = warpbreaks, corstr = "AR-M"
#' )
#' gee::qic(fit, type = "pan")
#' gee::qic(fit, type = "hh")
#' gee::qic(fit, type = "u")
#'
#' @export qic
qic <- function(
  object, type = "pan", tol = .Machine$double.eps, envir = parent.frame()
) {
  # Checking function parameter types ####
  if (!("gee" %in% class(object))) {
    stop("object must be an object produced by gee::gee")
  }
  if (!(type %in% c("pan", "hh", "u"))) {
    stop("type must be character value of the following options: pan, hh, u")
  }
  if (!(class(tol) %in% c("numeric", "integer"))) {
    stop("tol must be a numeric or an integer value")
  }
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

  # Calculating QIC ####

  ## Identifying the number of mean parameters ####
  nparms <- length(object$coefficients)

  ## Calculating penalty term ####
  penalty <- ifelse(
    test = type != "u",
    yes = gee::cic(object = object, type = type, tol = tol, envir = envir),
    no = nparms
  )

  ## Calculating QIC ####
  qic <- (-2 * gee::quasilik(object = object)) + (2 * penalty)

  # Return output
  return(qic)
}

#' @title Quasi-likelihood under the Indepence Model
#'
#' @description Calculating the quasi-likelihood under the independence model
#' based on estimating equations proposed by Liang and Zeger (1986).
#'
#' @param object An object of class \code{"gee"} representing the fit.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, 73 13-22.
#'
#' @return A numeric vector containing the quasi-likelihood under the
#' independence model value for the \code{"gee"} object under consideration.
#'
#' @examples
#' # Marginal analysis of random effects model for wool
#' data(warpbreaks)
#' fit <- gee::gee(
#'  breaks ~ tension, id = wool, data = warpbreaks, corstr = "AR-M"
#' )
#' gee::quasilik(object = fit)
#'
#' @export quasilik
quasilik <- function(object) {

  # Checking function parameter types ####
  if (!("gee" %in% class(object))) {
    stop("object must be an object produced by gee::gee")
  }

  # Obtaining family information ####
  fam <- object$family

  # Calculating quasi-likelihood ####
  y <- object$y
  mu <- object$fitted.values
  if (fam$family == "poisson") {
    quasi <- sum((y * log(mu)) - mu)
  } else if (fam$family == "gaussian") {
    quasi <- sum(((y - mu)^2) / -2)
  } else if (fam$family == "binomial") {
    quasi <- sum(y * log(mu / (1 - mu)) + log(1 - mu))
  } else if (fam$family == "Gamma") {
    quasi <- sum(-y / (mu - log(mu)))
  } else {
    stop("Error: distribution not recognized")
  }

  # Return output
  return(quasi)
}