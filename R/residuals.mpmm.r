##' Extract one-step-ahead residuals
##'
##' @title residuals
##' @param object an mpmm fit object
##' @param method character naming the method to calculate one-step-ahead
##'  residuals
##' @param trace logical; print progress to console
##' @param parallel logical; compute in parallel
##' @param ncores integer; number of cores to use (default = total cores
##' detected - 1)
##' @param ... additional arguments to be ignored
##'
##' @importFrom tibble tibble
##' @importFrom parallel detectCores
##' @method residuals mpmm
##' @return a list with components
##' \item{\code{res}}{a tibble with one-step-ahead residuals for longitude and
##'  latitude}
##'
##' @details Wrapper function for modified \code{\link{oneStepPredict}}
##' that calculates one-step-ahead residuals, which are residuals that account
##' for temporal correlation in latent states. The modification allows easier
##' parallel computation.
##'
##' @references Thygesen, U. H., C. M. Albertsen, C. W. Berg, K. Kristensen, and
##' A. Neilsen. 2017. Validation of ecological state space models using the
##' Laplace approximation. Environmental and Ecological Statistics 24:317–339.
##'
##' @examples
##' \dontrun{
##' data(ellie.ice)
##' fit <- mpmm(~ ice + (1 | id), data = ellie.ice)
##' summary(fit)
##' residuals(fit)
##' }
##' @export

residuals.mpmm <- function(object, method="oneStepGaussianOffMode", trace = FALSE, parallel = TRUE, ncores = detectCores() - 1, ...) {

    if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

  kidx <- as.vector(t(matrix(1:length(cbind(object$data$lon, object$data$lat)), ncol=2)))
  if(parallel)  sprintf("calculating residuals in parallel across %d cores...", ncores)
  mpmmres <- TMBoneStepPredict(object$tmb, observation.name ="ll",
                            data.term.indicator = "keep",
                            method = method,
                            discrete = FALSE,
                            subset = kidx,
                            trace = trace,
                            parallel = parallel,
                            ncores = ncores)
  res <-  tibble(id = object$data$id, date = object$data$date, res.lon = mpmmres[seq(1, nrow(mpmmres), by=2), "residual"],
                res.lat = mpmmres[seq(2, nrow(mpmmres), by=2), "residual"])

  class(res) <- append("mpmm_resid", class(res))

  return(res)

}
