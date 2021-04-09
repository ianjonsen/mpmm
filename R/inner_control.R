##' \code{inner_control} allows specification of control parameters for the
##' inner (newton) optimization used by \code{mpmm}
##'
##' See \code{\link{TMB::MakeADFun}} and \code{\link{TMB::newton}}
##' for details and available options. Adapted from S. Wotherspoon
##' \url{https://github.com/SWotherspoon/RWalc/blob/master/R/RWalc.R}
##' @title Control Values for \code{mpmm}
##' @param ... control parameters for the inner (newton) optimizer
##' @return Returns a list with components:
##'  \item{\code{control}}{list of control parameters for inner optimizer}
##' @seealso \code{\link{TMB::MakeADFun}}, \code{\link{TMB::newton}}
##' @examples
##' fit <- mpmm(~ ice + (ice | id),
##' data = ellie.ice.short,
##' inner_control = list(
##'     maxit = 1500,
##'     tol = 1e-04))
##' @export
inner_control <- function(...) {

  dots <- list(...)
  pars <- list(maxit = 1000,
               tol = 1e-08,
               smartsearch = TRUE)

  ## Override control parameters
  pars[names(dots)] <- dots
  list(in_control = pars)
}
