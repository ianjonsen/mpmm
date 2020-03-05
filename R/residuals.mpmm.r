##' Extract one-step-ahead residuals
##'
##' @title residuals
##' @param m a fitted object of class mpmm
##' @param method Character naming the method to calculate one-step-ahead residuals
##'
##' @importFrom TMB oneStepPredict
##' @importFrom tibble tibble
##' @method residuals mpmm
##' @return a list with components
##' \item{\code{res}}{a tibble with one-step-ahead residuals for longitude and latitude}
##'
##' @details Wrapper function for TMB::oneStepPredict function that calculates one-step-ahead residuals, which are residuals that account for temporal correlation in latent states.
##'
##' @references Thygesen, U. H., C. M. Albertsen, C. W. Berg, K. Kristensen, and A. Neilsen. 2017. Validation of ecological state space models using the Laplace approximation. Environmental and Ecological Statistics 24:317–339.
##'
##' @examples
##' \dontrun{
##' data(ellie.ice)
##' fit <- mpmm(~ ice + (1 | id), data = ellie.ice)
##' summary(fit)
##' residuals(fit)
##' }
##' @export

residuals.mpmm <- function(m, method="oneStepGaussianOffMode") {

  kidx <- as.vector(t(matrix(1:length(cbind(m$data$lon, m$data$lat)), ncol=2)))
  mpmmres <- oneStepPredict(m$tmb, observation.name ="ll",
                            data.term.indicator = "keep",
                            method = method,
                            discrete = FALSE,
                            subset = kidx, trace=FALSE)
  res <-  tibble(id = m$data$id, date = m$data$date, res.lon = mpmmres[seq(1, nrow(mpmmres), by=2), "residual"],
                res.lat = mpmmres[seq(2, nrow(mpmmres), by=2), "residual"])
  return(res)

}
