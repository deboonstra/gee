#' @title Correlation information criteria
#'
#' @description Calculating the correlation information criteria (CIC) based on
#' Hin and Wang (2009) or Hardin and Hilbe (2012).
#'
#' @param object An object of class \code{"gee"} representing the fit.
#' @param type A character string specifying the type of CIC to calculate. The
#' following are permitted: \code{"pan"} and \code{"hh"}.
#' Default is \code{"pan"}.
#' @param tol A numeric value specifying the tolerance when calculating the
#' inverse of the model-based covariance matrix.
#' Default is \code{.Machine$double.eps}.
#' @param envir An environment object to which find the \code{"gee"} object.
#'
#' @details The calculation of CIC proposed by Hin and Wang (2009) is a
#' modification of the quasi-likehood under the working independence model
#' (QIC) criteria, which was proposed by Pan (2001a), by only using the penalty
#' term of QIC. This version of CIC uses the mean structure parameter estimates
#' and the scale parameter estimate based on the working correlation structure
#' to calculate the model-based covariance matrix assuming the correlation
#' structure is independent. However, in 2012, Hardin and Hilbe proposed to
#' calculate CIC using the model-based covariance matrix given the mean
#' structure parameter estimates and the scale parameter estimate produced by
#' the independence model. Hin and Wang noted that the difference between the
#' two calculations of CIC is \eqn{O_{p}(n^{-1/2})}; thus, they do not expect
#' there to be much of a difference between the two variants of CIC.
#' To calculate the CIC value using the method proposed by Hin and Wang specify
#' \code{type = "pan"}, while to get the CIC value using the method proposed by
#' Hardin and Hilbe specify \code{type = "hh"}.
#'
#' If the MASS package is loaded then the \code{\link{ginv}} function is used
#' for matrix inversion. Otherwise the standard \code{\link{solve}} function is
#' used.
#'
#' @return A named numeric vector containing the CIC value for the \code{"gee"}
#' object under consideration. The name of the value is based on \code{type},
#' where \code{"CIC"} denotes \code{type = "pan"} and \code{"CIChh"} denotes
#' \code{type = "hh"}.
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
#' gee::cic(fit, type = "pan")
#' gee::cic(fit, type = "hh")
#'
#' @export cic
cic <- function(
  object, type = "pan", tol = .Machine$double.eps, envir = parent.frame()
) {
  # Checking function parameter types ####
  if (!("gee" %in% class(object))) {
    stop("object must be an object produced by gee::gee")
  }
  if (!(type %in% c("pan", "hh"))) {
    stop("type must be character value of the following options: pan, hh")
  }
  if (!(class(tol) %in% c("numeric", "integer"))) {
    stop("tol must be a numeric or an integer value")
  }
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

  # Setup functions
  invert <- if ("MASS" %in% loadedNamespaces()) {
    MASS::ginv
  } else {
    solve
  }

  # Obtaining covariance matrices ####

  ## Pulling modeling call to re-fit models ####
  model_call <- object$call

  ## Assigning corstr to independence ####
  model_call$corstr <- "independence"

  ## Re-fitting model based on CIC type ####
  if (type == "pan") {
    ### Forcing initial coefficient estimates to be final estimates of object
    model_call$b <- object$coefficients
    ### Forcing scale to be fixed based on working correlation structure
    model_call$scale.fix <- TRUE
    model_call$scale.value <- object$scale

    ### Re-fitting model
    model_indep <- suppressWarnings(eval(model_call, envir = envir))
  } else {
    ### Re-fitting model
    model_indep <- eval(model_call, envir = envir)
  }

  ## Getting covariance matrices ####

  ### Inverse of model-based covariance under independence ####
  omega <- invert(model_indep$naive.variance, tol = tol)

  ### Robust covariance given working correlation structure ####
  vr <- object$robust.variance

  # Calculating CIC ####
  trace <- sum(diag(x = crossprod(x = omega, y = vr)))
  if (type == "hh") {
    names(trace) <- "CIChh"
  } else {
    names(trace) <- "CIC"
  }

  # Returning CIC
  return(trace)
}